/**
 * Most of this file is straight out of the Next.JS Custom Server section -- it runs the nextjs app
 */
import { createServer } from 'http'
import { parse } from 'url'
import next from 'next'
import path from 'path'
import conf from '@/app/next.config'
import { io } from 'socket.io-client'
import * as dict from '@/utils/dict'

const dir = path.join(path.dirname(__dirname), 'app')
const dev = process.env.NODE_ENV !== 'production'
const hostname = '0.0.0.0'
const port = 3000

const [_node, _script, config] = process.argv
const {
  url,
  session_id,
  auth_token,
  project,
} = JSON.parse(config)
process.env.UFS_STORAGE = JSON.stringify({
  "cls": "ufs.impl.prefix.Prefix",
  "ufs": {
    "cls": "ufs.impl.sync.Sync",
    "ufs": {
      "cls": "ufs.impl.sbfs.SBFS",
      "auth_token": auth_token,
      "api_endpoint": "https://cavatica-api.sbgenomics.com",
      "drs_endpoint": "drs://cavatica-ga4gh-api.sbgenomics.com",
      "ttl": 60,
    },
  },
  "prefix": `/${project}`,
})
process.env.N_WORKERS = '50'
process.env.NEXTAUTH_SECRET = auth_token

const ctx = {
  lastMessage: Date.now(),
}
setInterval(() => {
  if ((Date.now() - ctx.lastMessage) > 2*60*1000) {
    console.log(`Session expired, exiting...`)
    process.exit(0)
  }
}, 30*1000)

console.log(`Connecting to ${url}...`)
fetch(`${url}/api/socket`).then(() => {
  const socket = io(url) // e.g. ws://localhost:3000
  socket.on('connect', () => {
    console.log(`Connected, joining ${session_id}...`)
    socket.emit('cavatica:join', session_id)
  })
  socket.on('http:send', async ({ id, path, headers, method, body }: { id: string, path: string, headers: Record<string, string>, method: string, body?: string }) => {
    console.log(JSON.stringify({ handle: { id, path, headers, method } }))
    ctx.lastMessage = Date.now()
    let responseHeaders: Record<string, string> = {}
    try {
      const req = await fetch(`http://${hostname}:${port}${path}`, { headers, method, body: body ? Buffer.from(body, 'base64') : undefined })
      responseHeaders = dict.fromHeaders(req.headers)
      const res = await req.text()
      const status = req.status
      socket.emit(`http:recv`, { id, status, body: res, headers: responseHeaders })
    } catch (err) {
      const status = 500
      console.warn(err)
      const res = JSON.stringify(err)
      socket.emit(`http:recv`, { id, status, body: res, headers: responseHeaders })
    }
  })
  socket.on('cavatica:close', () => {
    console.log(`Room has closed, exiting...`)
    socket.close()
    process.exit(0)
  })
})

// when using middleware `hostname` and `port` must be provided below
const app = next({ dev, hostname, port, dir, conf })
const handle = app.getRequestHandler()

app.prepare().then(() => {
  createServer(async (req, res) => {
    try {
      // Be sure to pass `true` as the second argument to `url.parse`.
      // This tells it to parse the query portion of the URL.
      if (!req.url) throw new Error('url is undefined')
      const parsedUrl = parse(req.url, true)
      await handle(req, res, parsedUrl)
    } catch (err) {
      console.error('Error occurred handling', req.url, err)
      res.statusCode = 500
      res.end('internal server error')
    }
  })
    .once('error', (err) => {
      console.error(err)
      process.exit(1)
    })
    .listen(port, () => {
      console.log(`> Ready on http://${hostname}:${port}`)
    })
}).catch((e) => console.error(e))
