import { parse } from 'url'
import http from 'http'
import next from 'next'
import path from 'path'
import { EventEmitter } from 'events'
import { isClaimedUpgrade } from '../upgrade'
import conf from '@/app/next.config'

export default async function plugin(server: http.Server, opts: { hostname: string, port: number }) {
  console.log('Starting NextJS server...')
  const dir = process.env.APP_ROOT ? `${process.env.APP_ROOT}/app` : path.dirname(path.dirname(__dirname))
  const dev = process.env.NODE_ENV !== 'production'
  // next.js attaches its upgrade handler to `httpServer`, give it a relay which only
  //  receives the upgrades other plugins (socket.io, tRPC) haven't claimed
  const nextUpgrades = new EventEmitter()
  server.on('upgrade', (req, socket, head) => {
    if (!isClaimedUpgrade(req)) nextUpgrades.emit('upgrade', req, socket, head)
  })
  const nextApp = next({ dev, hostname: opts.hostname, port: opts.port, dir, conf, httpServer: nextUpgrades as unknown as http.Server })
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
