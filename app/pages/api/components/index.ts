import handler from '@/utils/next-rest'
import krg from "@/app/krg"
import { metanodeType } from '@/app/fragments/components/listing'
import { UnsupportedMethodError } from '@/spec/error'

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const components = [
    ...krg.getDataNodes(),
    ...krg.getProcessNodes(),
  ].map((metanode) => {
    const { spec: id, meta: { icon: _, ...meta }, kind } = metanode
    return {
      id,
      kind,
      type: metanodeType(metanode),
      meta,
    }
  })
  res.status(200).json(components)
})
