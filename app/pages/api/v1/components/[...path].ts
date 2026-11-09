/**
 * Example files used to be served by component routes at these URLs, and those
 *  URLs were recorded as the source of uploads registered from the examples
 *  (`upload.url`, resolved when a DRS object is read). Keep serving them so that
 *  existing workflows built from example files keep working.
 *
 * New example uploads use the `components.file.register` procedure instead.
 */
import fs from 'fs'
import path from 'path'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'
import { legacyApiAssets } from '@/components/legacy-api'

export default handler(async (req, res) => {
  if (req.method !== 'GET' && req.method !== 'HEAD') throw new UnsupportedMethodError(req.method)
  const segments = Array.isArray(req.query.path) ? req.query.path : [req.query.path ?? '']
  const asset = legacyApiAssets[segments.filter(Boolean).join('/')]
  if (asset === undefined) throw new NotFoundError()
  const file = path.resolve(process.env.APP_ROOT as string, 'components', asset)
  const stat = await fs.promises.stat(file).catch(() => undefined)
  if (!stat?.isFile()) throw new NotFoundError()
  res.setHeader('Content-Length', stat.size)
  if (req.method === 'HEAD') return void res.status(200).end()
  await new Promise<void>((resolve, reject) => {
    fs.createReadStream(file).on('error', reject).pipe(res).on('finish', () => resolve())
  })
})
