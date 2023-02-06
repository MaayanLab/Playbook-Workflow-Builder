import { ResponseCodedError } from '@/spec/error'
import type { NextApiRequest, NextApiResponse } from 'next'

const handler = (handler_: (req: NextApiRequest, res: NextApiResponse) => Promise<void>) => async (req: NextApiRequest, res: NextApiResponse) => {
  try {
    await handler_(req, res)
  } catch (e) {
    res
      .status(('error_code' in (e as ResponseCodedError)) ? (e as ResponseCodedError).error_code : 500)
      .end((e as Error).toString())
  }
}
export default handler
