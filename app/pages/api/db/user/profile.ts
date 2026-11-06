import db from '@/app/db'
import { z } from 'zod'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import handler from '@/utils/next-rest'
import { NotFoundError, UnauthorizedError, UnsupportedMethodError } from '@/spec/error'

const BodyType = z.object({
  name: z.string().optional(),
  image: z.string().optional(),
  affiliation: z.string().optional(),
}).strict()

export default handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  if (req.method === 'GET') {
    const user = await db.selectFrom('user').selectAll().where('id', '=', session.user.id).executeTakeFirst()
    if (user === undefined) throw new NotFoundError()
    return res.status(200).json({
      image: user.image,
      name: user.name,
      email: user.email,
      affiliation: user.affiliation,
    })
  } else if (req.method === 'POST') {
    const data = BodyType.parse(req.body)
    const query = Object.keys(data).length > 0
      ? db.updateTable('user').set(data).where('id', '=', session.user.id).returningAll()
      : db.selectFrom('user').selectAll().where('id', '=', session.user.id)
    return res.status(200).json(await query.executeTakeFirst() ?? null)
  } else {
    throw new UnsupportedMethodError(req.method)
  }
})
