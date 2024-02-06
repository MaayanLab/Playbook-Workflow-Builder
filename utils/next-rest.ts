import { ResponseCodedError } from '@/spec/error'
import getRawBody from 'raw-body'
import type { NextApiRequest, NextApiResponse } from 'next'

export type NextApiRequestWithExtras = NextApiRequest & { json: () => Promise<any> }
export type RouteHandler = ReturnType<typeof handler>

function augmentRequestWithExtras(req: NextApiRequest) {
  async function json() {
    const body = await getRawBody(req)
    return JSON.parse(body.toString())
  }
  Object.assign(req, { json })
  return req as NextApiRequestWithExtras
}

const handler = (handler_: (req: NextApiRequestWithExtras, res: NextApiResponse) => Promise<void>) => async (req: NextApiRequest, res: NextApiResponse) => {
  try {
    await handler_(augmentRequestWithExtras(req), res)
  } catch (e) {
    res
      .status(('error_code' in (e as ResponseCodedError)) ? (e as ResponseCodedError).error_code : 500)
      .send((e as ResponseCodedError).message ?? (e as Error).toString())
  }
}
export default handler
