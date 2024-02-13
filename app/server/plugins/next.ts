import http from 'http';
import { parse } from 'url'
import next from 'next'
import path from 'path'
import conf from '@/app/next.config'

export default async function plugin(opts: { hostname: string, port: number }) {
  console.log('Starting NextJS server...')
  const dir = path.dirname(path.dirname(__dirname))
  const dev = process.env.NODE_ENV !== 'production'
  const app = next({ dev, hostname: opts.hostname, port: opts.port, dir, conf })
  const handle = app.getRequestHandler()
  await app.prepare()
  const server = http.createServer(async (req, res) => {
    try {
      // Be sure to pass `true` as the second argument to `url.parse`.
      // This tells it to parse the query portion of the URL.
      if (!req.url) throw new Error('url is undefined')
      const parsedUrl = parse(req.url, true)
      if (parsedUrl.pathname?.startsWith('/socket.io')) {
        console.log('ignoreing socket.io')
        return
      }
      await handle(req, res, parsedUrl)
    } catch (err) {
      console.error('Error occurred handling', req.url, err)
      res.statusCode = 500
      res.end('internal server error')
    }
  })
  server.once('error', (err) => {
    console.error(err)
    process.exit(1)
  })
  await new Promise<void>((resolve, reject) => {
    server.listen(opts.port, () => {resolve()})
  })
  console.log(`Server listening on http://${opts.hostname}:${opts.port}`);
}
