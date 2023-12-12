import { UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'
import * as dict from '@/utils/dict'

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  res.status(200).json(
    dict.filter(
      process.env as Record<string, string>,
      ({ key }) => key.startsWith('NEXT_PUBLIC_'),
    )
  )
})