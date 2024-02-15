import { parse } from 'url'
import http from 'http'
import next from 'next'
import path from 'path'
import conf from '@/app/next.config'

export default async function plugin(server: http.Server, opts: { hostname: string, port: number }) {
  console.log('Starting NextJS server...')
  const dir = process.env.APP_ROOT ? `${process.env.APP_ROOT}/app` : path.dirname(path.dirname(__dirname))
  const dev = process.env.NODE_ENV !== 'production'
  const nextApp = next({ dev, hostname: opts.hostname, port: opts.port, dir, conf })
  const nextHandler = nextApp.getRequestHandler()
  await nextApp.prepare()
  server.on('request', async (req, res) => {
    try {
      if (!req.url) throw new Error('Unhandled Error')
      const parsedUrl = parse(req.url, true)
      await nextHandler(req, res, parsedUrl)
    } catch (err) {
      console.error('Error occurred handling', req.url, err)
      res.statusCode = 500
      res.end('internal server error')
    }
  })
}
