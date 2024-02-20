import { z } from 'zod'
import http from 'http'

export const Options = z.object({
  hostname: z.string().default('0.0.0.0'),
  port: z.number().default(3000),
  plugins: z.array(z.enum(['ws', 'next', 'cavatica-proxy'])).default(['ws', 'next']).transform(L => new Set(L)),
  proxy: z.object({
    url: z.string(),
    session_id: z.string(),
    auth_token: z.string(),
    project: z.string(),
  }).optional(),
})

export default async function main(opts: z.infer<typeof Options>) {
  const server = http.createServer()
  if (opts.plugins.has('ws')) {
    await (await import('./plugins/ws')).default(server, opts)
  }
  if (opts.plugins.has('next')) {
    await (await import('./plugins/next')).default(server, opts)
  }
  if (opts.plugins.has('cavatica-proxy')) {
    await (await import('./plugins/cavatica-proxy')).default(server, opts)
  }
  server.once('error', (err) => {
    console.error(err)
    process.exit(1)
  })
  await new Promise<void>((resolve, reject) => {
    server.listen(opts.port, () => {resolve()})
  })
  console.log(`Server listening on http://${opts.hostname}:${opts.port}`);
}
