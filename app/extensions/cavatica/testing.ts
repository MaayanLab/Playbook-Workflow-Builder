/**
 * A stand-in for @/app/extensions/cavatica which runs the worker on this machine
 *  instead of in a CAVATICA task, enabled with CAVATICA_LOCAL_WORKER=1.
 */
import { spawn } from "child_process"

async function sleep(s: number) {
  return await new Promise<void>((resolve, reject) => {setTimeout(() => resolve(), s*1000)})
}

export async function *run_wes_worker({
  url,
  session_id,
  auth_token,
  project,
  api_endpoint = 'https://cavatica-api.sbgenomics.com',
  wes_endpoint = 'wes://cavatica-ga4gh-api.sbgenomics.com',
  polling_interval = 5,
}: {
  url: string,
  session_id: string,
  auth_token: string,
  project: string,
  api_endpoint?: string,
  wes_endpoint?: string,
  polling_interval?: number,
}) {
  // run the worker locally (with its own in-memory database & `next dev` build directory)
  // DATABASE_URL is set (empty) rather than removed, next.config.js would otherwise load it from .env
  const proc = spawn('npx', [
    'tsx', 'cli/cavatica-worker.ts',
    JSON.stringify({ url, session_id, auth_token, project }),
  ], {
    cwd: process.env.APP_ROOT,
    env: { ...process.env, DATABASE_URL: '', PORT: '3001', CAVATICA_WORKER_NEXT: 'dev', NEXT_DIST_DIR: '.next-cavatica-worker' },
  })
  const ctx = {
    current_state: null as string | null,
    state: null as string | null,
  }
  proc.stdout.on('data', (msg) => console.log(`worker: ${Buffer.from(msg).toString()}`))
  proc.stderr.on('data', (msg) => console.error(`worker: ${Buffer.from(msg).toString()}`))
  proc.on('close', () => {ctx.current_state = 'COMPLETE'})
  // like the WES task, report the run id once started (it's needed to disconnect)
  yield { run_id: '12345', state: null }
  while (true) {
    await sleep(polling_interval * (0.5 + Math.random()))
    if (ctx.state !== ctx.current_state) {
      ctx.state = ctx.current_state
      yield { run_id: '12345', state: ctx.state }
    }
    if (ctx.state === 'COMPLETE') {
      break
    }
  }
}

export async function abort_wes_worker({
  run_id,
  auth_token,
  wes_endpoint = 'wes://cavatica-ga4gh-api.sbgenomics.com',
}: { run_id: string, auth_token: string, wes_endpoint?: string }) {
  return true
}
