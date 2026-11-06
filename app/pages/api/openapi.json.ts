import type { NextApiRequest, NextApiResponse } from 'next'
import { openApiDocument } from '@/app/openapi'

export default function handler(_req: NextApiRequest, res: NextApiResponse) {
  res.status(200).json(openApiDocument())
}
