import db from '@/app/db'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import handler from "@/utils/next-rest"

export default handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  if (req.method === 'GET') {
    const suggestions = await db.objects.suggestion.findMany({ where: { user: session.user.id } })
    return res.status(200).json(suggestions)
  } else {
    throw new UnsupportedMethodError(req.method)
  }
})
