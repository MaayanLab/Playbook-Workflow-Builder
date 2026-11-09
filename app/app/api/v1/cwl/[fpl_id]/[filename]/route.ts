import YAML from 'yaml'
import fpprg from '@/app/fpprg'
import { getKRG } from '@/app/server-krg'
import { cwl_for_playbook } from '@/core/fpl2cwl'
import { NotFoundError } from '@/spec/error'
import { routeResponse } from '@/app/app/api/_lib/response'

export async function GET(
  _request: Request,
  props: { params: Promise<{ fpl_id: string, filename: string }> }
) {
  const params = await props.params;
  return routeResponse(async () => {
    const fpl = await fpprg.getFPL(params.fpl_id)
    if (!fpl) throw new NotFoundError()
    const files = await cwl_for_playbook({ krg: await getKRG(), fpl })
    const file = files[params.filename as keyof typeof files]
    if (!file) throw new NotFoundError()
    return new Response(YAML.stringify(file), { headers: { 'Content-Type': 'application/yaml' } })
  })
}
