import { ResponseCodedError } from '@/spec/error'
import type { NextApiRequest, NextApiResponse } from 'next'

const handler = (handler_: (req: NextApiRequest, res: NextApiResponse) => Promise<void>) => async (req: NextApiRequest, res: NextApiResponse) => {
  try {
    await handler_(req, res)
  } catch (e) {
    res.status(
      (e instanceof ResponseCodedError) ? e.error_code : 500
    ).end((e as Error).toString())
  }
}
export default handler
