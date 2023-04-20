import handler from '@/utils/next-rest'
import db from '@/app/db'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { NotFoundError, UnauthorizedError } from '@/spec/error'
import { z } from 'zod'

const QueryType = z.object({
  object_id: z.string(),
  access_id: z.string(),
})
const BodyType = z.object({
  pasports: z.array(z.string()).optional(),
})

export default handler(async (req, res) => {
  const { object_id } = QueryType.parse(req.query)
  if (req.method === 'POST') {
    const body = BodyType.parse(req.body)
    // TODO: potentially authorize with passports
  }
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  const drs_access = await db.objects.drs_access.findUnique({ where: { access_id, object_id } })
  if (drs_access === null) throw new NotFoundError()
  res.json({
    "url": drs_access.url,
    "headers": drs_access.headers,
  })
})
