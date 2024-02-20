import { UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  res.status(200).json({
    uri: process.env.PUBLIC_URL,
    transports: ['websocket'],
  })
})
