import fpprg from '@/app/fpprg'
import { getKRG } from '@/app/server-krg'
import { generateExport } from '@/app/exports'
import { NotFoundError } from '@/spec/error'
import { openApiFetchHandler } from '@/app/openapi'
import { attachmentHeaders, routeResponse } from '@/app/app/api/_lib/response'

export async function GET(request: Request, props: { params: Promise<{ fpl_id: string }> }) {
  const params = await props.params
  return routeResponse(async () => {
    if ((new URL(request.url).searchParams.get('format') ?? 'zip') === 'json') return openApiFetchHandler(request)
    const fpl = await fpprg.getFPL(params.fpl_id)
    if (!fpl) throw new NotFoundError()
    const file = await generateExport({ krg: await getKRG(), fpl, format: 'ro-crate-zip', params: {} })
    return new Response(file.data, { headers: attachmentHeaders(file.filename, file.data.length) })
  })
}
