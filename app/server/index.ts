import { z } from 'zod'

export const Options = z.object({
  next: z.object({
    hostname: z.string().default('0.0.0.0'),
    port: z.number().default(3000),
  }).optional(),
  ws: z.object({
    hostname: z.string().default('0.0.0.0'),
    port: z.number().default(3005),
  }).optional(),
})

export default async function main(opts: z.infer<typeof Options>) {
  if (opts.ws) {
    await (await import('./plugins/ws')).default(opts.ws)
  }
  if (opts.next) {
    await (await import('./plugins/next')).default(opts.next)
  }
}
