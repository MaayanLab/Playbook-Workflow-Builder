import fs from 'fs'
import path from 'path'
import { legacyApiAssets } from '@/components/legacy-api'
import { NotFoundError } from '@/spec/error'
import { nodeStreamResponse, routeResponse } from '@/app/app/api/_lib/response'

async function asset(params: { path: string[] }) {
  const relative = legacyApiAssets[params.path.join('/')]
  if (!relative) throw new NotFoundError()
  const file = path.resolve(process.env.APP_ROOT as string, 'components', relative)
  const stat = await fs.promises.stat(file).catch(() => undefined)
  if (!stat?.isFile()) throw new NotFoundError()
  return { file, stat }
}

export async function HEAD(_request: Request, props: { params: Promise<{ path: string[] }> }) {
  const params = await props.params;
  return routeResponse(async () => {
    const { stat } = await asset(params)
    return new Response(null, { headers: { 'Content-Length': String(stat.size) } })
  })
}

export async function GET(_request: Request, props: { params: Promise<{ path: string[] }> }) {
  const params = await props.params;
  return routeResponse(async () => {
    const { file, stat } = await asset(params)
    return nodeStreamResponse(fs.createReadStream(file), { headers: { 'Content-Length': String(stat.size) } })
  })
}
