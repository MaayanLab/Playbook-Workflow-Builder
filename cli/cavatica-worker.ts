/**
 * CAVATICA worker: runs the app inside a CAVATICA task and serves a session's
 *  requests through the public server's tunnel (see @/app/extensions/cavatica/tunnel).
 *
 * Usage: tsx cli/cavatica-worker.ts '{"url":"https://public.server","session_id":"...","auth_token":"...","project":"..."}'
 *
 * Environment:
 *  - PORT: the local port to run the app on (default 3001)
 *  - CAVATICA_WORKER_NEXT=dev: run the app with `next dev` instead of `next start` (local testing,
 *     use with NEXT_DIST_DIR to keep its build separate, see @/app/extensions/cavatica/testing)
 *  - CAVATICA_WORKER_APP_URL: use an already running app instead of starting one (for testing)
 */
import { spawn, type ChildProcess } from 'child_process'
import path from 'path'
import { z } from 'zod'
import { EventSource } from 'eventsource'
import { createTRPCClient, httpLink, httpSubscriptionLink, splitLink } from '@trpc/client'
import type { AppRouter } from '@/app/trpc/router'
import type { WorkerEvent, WorkerResponse } from '@/app/extensions/cavatica/tunnel'

const ConfigC = z.object({
  url: z.string(),
  session_id: z.string(),
  auth_token: z.string(),
  project: z.string(),
})

/** exit after this long without requests (and none in flight) */
const IDLE_TIMEOUT_MS = 2 * 60 * 1000
/** responses are sent in chunks of up to this size */
const CHUNK_SIZE = 256 * 1024
/** request headers which only apply to the connection to the public server */
const hopByHop = new Set(['connection', 'keep-alive', 'transfer-encoding', 'upgrade', 'host', 'content-length', 'accept-encoding'])

const config = ConfigC.parse(JSON.parse(process.argv[2] ?? '{}'))
const publicUrl = config.url.replace(/\/$/, '')
const port = +(process.env.PORT || 3001)
const appUrl = process.env.CAVATICA_WORKER_APP_URL ?? `http://127.0.0.1:${port}`

/** Start the app (`next start`) configured to use the user's CAVATICA project */
function startApp(): ChildProcess {
  const appRoot = process.env.APP_ROOT ?? path.dirname(__dirname)
  const dev = process.env.CAVATICA_WORKER_NEXT === 'dev'
  return spawn('npx', ['next', dev ? 'dev' : 'start', '-p', `${port}`, '-H', '127.0.0.1'], {
    cwd: path.join(appRoot, 'app'),
    stdio: 'inherit',
    env: {
      ...process.env,
      APP_ROOT: appRoot,
      NODE_ENV: dev ? 'development' : 'production',
      // the worker uses its own in-memory database (set empty, next.config.js would load a missing one from .env)
      DATABASE_URL: '',
      UFS_STORAGE: JSON.stringify({
        "cls": "ufs.impl.prefix.Prefix",
        "ufs": {
          "cls": "ufs.impl.sync.Sync",
          "ufs": {
            "cls": "ufs.impl.sbfs.SBFS",
            "auth_token": config.auth_token,
            "api_endpoint": "https://cavatica-api.sbgenomics.com",
            "drs_endpoint": "drs://cavatica-ga4gh-api.sbgenomics.com",
            "ttl": 60,
          },
        },
        "prefix": `/${config.project}`,
      }),
      N_WORKERS: '50',
      // requests forwarded by the public server carry this as `Authorization: Token ...`
      NEXTAUTH_SECRET: config.auth_token,
      PUBLIC_URL: appUrl,
      NEXT_PUBLIC_URL: appUrl,
    },
  })
}

async function waitForApp() {
  for (let attempt = 0; attempt < 300; attempt++) {
    try {
      await fetch(`${appUrl}/api/trpc/system.runtimeConfig`)
      return
    } catch {
      await new Promise(resolve => setTimeout(resolve, 1000))
    }
  }
  throw new Error('App did not start')
}

const authorization = `Token ${config.auth_token}`
const trpc = createTRPCClient<AppRouter>({
  links: [
    splitLink({
      condition: (op) => op.type === 'subscription',
      true: httpSubscriptionLink({
        url: `${publicUrl}/api/trpc`,
        EventSource,
        eventSourceOptions: {
          fetch: (input, init) => fetch(input, { ...init, headers: { ...init?.headers, authorization } }),
        },
      }),
      false: httpLink({ url: `${publicUrl}/api/trpc`, headers: { authorization } }),
    }),
  ],
})

const inFlight = new Map<string, AbortController>()
let lastActivity = Date.now()

async function send(response: WorkerResponse) {
  await trpc.tunnel.respond.mutate({ session_id: config.session_id, response })
}

/** Perform a request against the app and stream the response back through the tunnel */
async function handle(event: Extract<WorkerEvent, { type: 'request' }>) {
  const controller = new AbortController()
  inFlight.set(event.id, controller)
  try {
    const res = await fetch(`${appUrl}${event.path}`, {
      method: event.method,
      headers: Object.fromEntries(Object.entries(event.headers).filter(([key]) => !hopByHop.has(key.toLowerCase()))),
      body: event.body !== undefined ? Buffer.from(event.body, 'base64') : undefined,
      signal: controller.signal,
      redirect: 'manual',
    })
    await send({ type: 'head', id: event.id, status: res.status, headers: Object.fromEntries(res.headers.entries()) })
    if (res.body) {
      // chunks are sent in order (one at a time)
      for await (const chunk of res.body as unknown as AsyncIterable<Uint8Array>) {
        for (let offset = 0; offset < chunk.length; offset += CHUNK_SIZE) {
          const data = Buffer.from(chunk.subarray(offset, offset + CHUNK_SIZE)).toString('base64')
          await send({ type: 'chunk', id: event.id, data })
        }
      }
    }
    await send({ type: 'end', id: event.id })
  } catch (error) {
    if (!controller.signal.aborted) {
      console.warn(`request ${event.method} ${event.path} failed`, error)
      await send({ type: 'error', id: event.id, message: String(error) }).catch(() => {})
    }
  } finally {
    inFlight.delete(event.id)
    lastActivity = Date.now()
  }
}

async function main() {
  const app = process.env.CAVATICA_WORKER_APP_URL ? undefined : startApp()
  app?.on('exit', (code) => {
    console.log(`App exited (${code}), exiting...`)
    process.exit(code ?? 1)
  })
  const shutdown = (reason: string) => {
    console.log(`${reason}, exiting...`)
    app?.kill()
    process.exit(0)
  }
  await waitForApp()
  console.log(`Connecting to ${publicUrl} for session ${config.session_id}...`)
  trpc.tunnel.serve.subscribe({ session_id: config.session_id }, {
    onStarted: () => { console.log('Connected') },
    onData: (event) => {
      lastActivity = Date.now()
      if (event.type === 'request') handle(event)
      else if (event.type === 'abort') inFlight.get(event.id)?.abort()
      else if (event.type === 'close') shutdown('Session closed')
    },
    onError: (error) => { console.error('Tunnel error', error) },
  })
  setInterval(() => {
    if (inFlight.size === 0 && Date.now() - lastActivity > IDLE_TIMEOUT_MS) shutdown('Session expired')
  }, 30 * 1000)
}

main().catch((error) => {
  console.error(error)
  process.exit(1)
})
