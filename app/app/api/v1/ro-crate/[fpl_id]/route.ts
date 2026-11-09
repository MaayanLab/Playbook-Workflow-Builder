import JSZip from 'jszip'
import fpprg from '@/app/fpprg'
import { getKRG } from '@/app/server-krg'
import { fpl2ro_crate } from '@/core/fpl2ro-crate'
import { NotFoundError } from '@/spec/error'
import { openApiFetchHandler } from '@/app/openapi'
import { routeResponse } from '@/app/app/api/_lib/response'

export async function GET(request: Request, props: { params: Promise<{ fpl_id: string }> }) {
  const params = await props.params;
  return routeResponse(async () => {
    if ((new URL(request.url).searchParams.get('format') ?? 'zip') === 'json') return openApiFetchHandler(request)
    const fpl = await fpprg.getFPL(params.fpl_id)
    if (!fpl) throw new NotFoundError()
    const files = await fpl2ro_crate({ krg: await getKRG(), fpl })
    const zip = new JSZip()
    for (const [name, value] of Object.entries(files)) zip.file(name, value)
    return new Response(await zip.generateAsync({ type: 'arraybuffer' }), {
      headers: { 'Content-Type': 'application/zip' },
    })
  })
}
