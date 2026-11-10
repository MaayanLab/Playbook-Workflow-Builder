/**
 * Produce a LaTeX bundle (`format=zip`) or PDF (`format=pdf`) of a playbook in the request.
 *  This can take a long time, the UI instead uses `exports.create` (see app/exports.ts).
 */
import fpprg from '@/app/fpprg'
import { getKRG } from '@/app/server-krg'
import { ExportMetadataC, exportAuthor, generateExport } from '@/app/exports'
import { getRequestSessionWithId } from '@/app/extensions/auth/server'
import { NotFoundError } from '@/spec/error'
import { attachmentHeaders, routeResponse } from '@/app/app/api/_lib/response'

export async function GET(request: Request, props: { params: Promise<{ fpl_id: string }> }) {
  const params = await props.params
  return routeResponse(async () => {
    const query = new URL(request.url).searchParams
    const format = query.get('format') === 'pdf' ? 'tex-pdf' : 'tex-zip'
    const metadataRaw = query.get('metadata')
    const metadata = metadataRaw ? ExportMetadataC.parse(JSON.parse(metadataRaw)) : undefined
    const fpl = await fpprg.getFPL(params.fpl_id)
    if (!fpl) throw new NotFoundError()
    const session = await getRequestSessionWithId(request)
    const author = await exportAuthor(session?.user?.id)
    const file = await generateExport({ krg: await getKRG(), fpl, format, params: { metadata, author } })
    return new Response(file.data, { headers: attachmentHeaders(file.filename, file.data.length) })
  })
}
