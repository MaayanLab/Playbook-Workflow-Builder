import handler from '@/utils/next-rest'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError, UnsupportedMethodError } from "@/spec/error"
import { z } from 'zod'

export const GET = handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  if (!process.env.ELYSIUM_USERNAME || !process.env.ELYSIUM_PASSWORD) throw new Error('Not Configured')
  const { uid } = z.object({ uid: z.string() }).parse(req.query)
  
  const charonSearchParams = new URLSearchParams()
  charonSearchParams.append('username', process.env.ELYSIUM_USERNAME)
  charonSearchParams.append('password', process.env.ELYSIUM_PASSWORD)
  const charonReq = await fetch(`https://maayanlab.cloud/charon/signpolicy?${charonSearchParams.toString()}`)
  const charonRes = await charonReq.json()

  const formData: Record<string, string> = {}
  formData['key'] = charonRes['uid'] + `/${uid}-\${filename}`
  formData['AWSAccessKeyId'] = charonRes['cid']
  formData['acl'] = 'private'
  formData['success_action_redirect'] = 'success.html'
  formData['policy'] = charonRes['policy']
  formData['signature'] = charonRes['signature']
  formData['Content-Type'] = 'application/octet-stream'
  res.json(formData)
})
