import db from '@/app/db'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import handler from '@/utils/next-rest'
import { NotFoundError, UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import { z } from 'zod'

const QueryType = z.object({
  id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'POST') throw new UnsupportedMethodError()
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  const { id } = QueryType.parse(req.query)
  const deletedUserUpload = await db.objects.user_upload.delete({ where: { id, user: session.user.id } })
  if (deletedUserUpload === null) throw new NotFoundError()
  res.status(200).json(deletedUserUpload)
})
