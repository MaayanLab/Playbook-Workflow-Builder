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
    const user = await db.objects.user.findUnique({ where: { id: session.user.id } })
    if (user === null) throw new NotFoundError()
    return res.status(200).json({
      image: user.image,
      name: user.name,
      email: user.email,
      affiliation: user.affiliation,
    })
  } else if (req.method === 'POST') {
    return res.status(200).json(await db.objects.user.update({
      where: { id: session.user.id },
      data: BodyType.parse(req.body),
    }))
  } else {
    throw new UnsupportedMethodError(req.method)
  }
})
