import handler from '@/utils/next-rest'
import db from '@/app/db'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { NotFoundError, UnauthorizedError } from '@/spec/error'
import { z } from 'zod'

const QueryType = z.object({
  object_id: z.string(),
})
const BodyType = z.object({
  expand: z.boolean().optional(),
  passports: z.array(z.string()).optional(),
})

export default handler(async (req, res) => {
  const { object_id } = QueryType.parse(req.query)
  if (req.method === 'POST') {
    const body = BodyType.parse(req.body)
    // TODO: potentially authorize with passports
  }
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  const drs_object = await db.objects.drs.findUnique({ where: { id: object_id } })
  if (drs_object === null) throw new NotFoundError()
  if (req.method === 'OPTIONS') {
    res.status(200)
    return
  }
  res.json({
    "id": drs_object.drs_id,
    "name": drs_object.name,
    "self_uri": `${(process.env.PUBLIC_URL||'').replace(/^https/, 'drs').replace(/\/$/, '')}/${drs_object.drs_id}`,
    "size": drs_object.size,
    "created_time": drs_object.created_time,
    "checksums": drs_object.checksums,
    // optional
    "updated_time": drs_object.updated_time,
    "version": drs_object.version,
    "mime_type": drs_object.mime_type,
    "access_methods": drs_object.access_methods,
    "description": drs_object.description,
    "aliases": drs_object.aliases,
  })
})
