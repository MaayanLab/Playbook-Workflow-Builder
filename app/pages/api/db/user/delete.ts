import db from '@/app/db'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import handler from '@/utils/next-rest'
import { NotFoundError, UnauthorizedError, UnsupportedMethodError } from '@/spec/error'

export default handler(async (req, res) => {
  if (req.method !== 'POST') throw new UnsupportedMethodError()
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  const deletedUser = await db.objects.user.delete({ where: { id: session.user.id } })
  if (deletedUser === null) throw new NotFoundError()
  res.status(200).json(deletedUser)
})
