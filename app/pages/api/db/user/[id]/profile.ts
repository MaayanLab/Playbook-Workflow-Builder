import db from '@/app/db'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import { z } from 'zod'

const QueryType = z.object({
  id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { id } = QueryType.parse(req.query)
  const user = await db.selectFrom('user').selectAll().where('id', '=', id).executeTakeFirst()
  if (user === undefined) throw new NotFoundError()
  res.status(200).json({
    image: user.image,
    name: user.name,
    email: user.email,
    affiliation: user.affiliation,
  })
})
