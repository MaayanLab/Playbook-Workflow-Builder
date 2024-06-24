import handler from '@/utils/next-rest'
import { metanodeType } from '@/app/fragments/components/listing'
import { z } from 'zod'
import { UnsupportedMethodError } from '@/spec/error'
import krg from '@/app/krg'

const QueryType = z.object({
  spec: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { spec } = QueryType.parse(req.query)
  return res.status(200).json(
    krg.getPrevProcess(spec).map(metanode => {
      const { spec: id, meta: { icon: _, ...meta }, kind } = metanode
      return {
        id,
        kind,
        type: metanodeType(metanode),
        meta,
      }
    })
  )
})
