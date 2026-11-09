import fpprg from '@/app/fpprg'
import { openApiFetchHandler } from '@/app/openapi'
import { routeResponse } from '@/app/app/api/_lib/response'
import { NotFoundError } from '@/spec/error'

export function GET(request: Request) {
  return routeResponse(() => openApiFetchHandler(request))
}

export async function HEAD(_request: Request, props: { params: Promise<{ fpl_id: string }> }) {
  const params = await props.params;
  return routeResponse(async () => {
    if (!(await fpprg.getFPL(params.fpl_id))) throw new NotFoundError()
    return new Response(null, { status: 200 })
  });
}
