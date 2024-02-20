import http from 'http'
import { io } from 'socket.io-client'
import * as dict from '@/utils/dict'
import { z } from 'zod'
import type { Options } from '..'

async function wsConnect(url: string) {
  console.log(`Fetching ws config from ${url}...`)
  const req = await fetch(`${url}/api/socket`)
  const { uri, ...opts } = await req.json()
  console.log(`> Connecting to ${uri}...`)
  return io(uri, opts)
}

export default async function plugin(server: http.Server, opts: z.infer<typeof Options>) {
  if (!opts.proxy) throw new Error('Missing proxy config')
  console.log('Starting CAVATICA proxy...')

  process.env.UFS_STORAGE = JSON.stringify({
    "cls": "ufs.impl.prefix.Prefix",
    "ufs": {
      "cls": "ufs.impl.sync.Sync",
      "ufs": {
        "cls": "ufs.impl.sbfs.SBFS",
        "auth_token": opts.proxy.auth_token,
        "api_endpoint": "https://cavatica-api.sbgenomics.com",
        "drs_endpoint": "drs://cavatica-ga4gh-api.sbgenomics.com",
        "ttl": 60,
      },
    },
    "prefix": `/${opts.proxy.project}`,
  })
  process.env.N_WORKERS = '50'
  process.env.NEXTAUTH_SECRET = opts.proxy.auth_token
  process.env.PUBLIC_URL = process.env.NEXT_PUBLIC_URL = `http://${opts.hostname}:${opts.port}`

  // monitor socket messages
  //  -- close if around 2 minutes elapsed with no messages
  const ctx = {
    lastMessage: Date.now(),
  }
  setInterval(() => {
    if ((Date.now() - ctx.lastMessage) > 2*60*1000) {
      console.log(`Session expired, exiting...`)
      process.exit(0)
    }
  }, 30*1000)

  const publicServerSocket = await wsConnect(opts.proxy.url)
  publicServerSocket.on('connect', () => {
    if (!opts.proxy) return
    console.log(`Connected, joining ${opts.proxy.session_id}...`)
    publicServerSocket.emit('worker:join', opts.proxy.session_id)
  })
  publicServerSocket.on('http:send', async ({ id, path, headers, method, body }: { id: string, path: string, headers: Record<string, string>, method: string, body?: string }) => {
    console.log(JSON.stringify({ handle: { id, path, headers, method } }))
    ctx.lastMessage = Date.now()
    let responseHeaders: Record<string, string> = {}
    try {
      const req = await fetch(`http://${opts.hostname}:${opts.port}${path}`, { headers, method, body: body ? Buffer.from(body, 'base64') : undefined })
      responseHeaders = dict.fromHeaders(req.headers)
      const res = await req.text()
      const status = req.status
      publicServerSocket.emit(`http:recv`, { id, status, body: res, headers: responseHeaders })
    } catch (err) {
      const status = 500
      console.warn(err)
      const res = JSON.stringify(err)
      publicServerSocket.emit(`http:recv`, { id, status, body: res, headers: responseHeaders })
    }
  })
  publicServerSocket.on('cavatica:close', () => {
    console.log(`Room has closed, exiting...`)
    publicServerSocket.close()
    process.exit(0)
  })

  server.on('listening', async () => {
    // create a bidirectional channel between the two servers over websocket
    //  in the public server, this is prefixed by `ws:{id}:` but locally
    //  this prefix is omitted
    const privateServerSocket = await wsConnect(`http://${opts.hostname}:${opts.port}`)
    publicServerSocket.onAny((evt, ...args) => {
      if (!opts.proxy) return
      if (evt.startsWith(`ws:`)) {
        console.log(`forwarding ${evt.slice(`ws:`.length)} from public to private`)
        ctx.lastMessage = Date.now()
        privateServerSocket.emit(evt.slice(`ws:`.length), ...args)
      }
    })
    privateServerSocket.onAny((evt, ...args) => {
      if (!opts.proxy) return
      console.log(`forwarding ${evt} from private to public`)
      publicServerSocket.emit(`ws:${evt}`, ...args)
    })
  })
}
