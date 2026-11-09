import JSZip from 'jszip'
import YAML from 'yaml'
import fpprg from '@/app/fpprg'
import { getKRG } from '@/app/server-krg'
import { cwl_for_playbook } from '@/core/fpl2cwl'
import { NotFoundError } from '@/spec/error'
import { openApiFetchHandler } from '@/app/openapi'
import { routeResponse } from '@/app/app/api/_lib/response'

export async function GET(request: Request, props: { params: Promise<{ fpl_id: string }> }) {
  const params = await props.params;
  return routeResponse(async () => {
    const format = new URL(request.url).searchParams.get('format') ?? 'zip'
    if (format === 'json' || format === 'list') return openApiFetchHandler(request)
    const fpl = await fpprg.getFPL(params.fpl_id)
    if (!fpl) throw new NotFoundError()
    const files = await cwl_for_playbook({ krg: await getKRG(), fpl })
    const zip = new JSZip()
    for (const [name, value] of Object.entries(files)) zip.file(name, YAML.stringify(value))
    return new Response(await zip.generateAsync({ type: 'arraybuffer' }), {
      headers: { 'Content-Type': 'application/zip' },
    })
  })
}
