/**
 * Most of this file is straight out of the Next.JS Custom Server section -- it runs the nextjs app
 */
import { createServer } from 'http'
import { parse } from 'url'
import next from 'next'
import path from 'path'
import { io } from 'socket.io-client'
import conf from '@/app/next.config'

const dir = path.join(path.dirname(__dirname), 'app')
const dev = process.env.NODE_ENV !== 'production'
const hostname = '0.0.0.0'
const port = 3001
// when using middleware `hostname` and `port` must be provided below
const app = next({ dev, hostname, port, dir, conf })
const handle = app.getRequestHandler()

app.prepare().then(() => {
  /* BEGIN client proxy */
  console.log('starting')
  const socket = io(process.argv[2])
  
  socket.on('http:send', ({ id, path, method, body }: { id: string, path: string, method: string, body: string }) => {
    fetch(path, { method, body })
      .then(res => socket.emit(`http:recv:${id}`, { status: res.status, body: res.body }))
      .catch(err => socket.emit(`http:recv:${id}`, { status: 500, body: err }))
  })
  /* END client proxy */
  
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
