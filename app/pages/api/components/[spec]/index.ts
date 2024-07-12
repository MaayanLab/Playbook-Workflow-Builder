import handler from '@/utils/next-rest'
import { metanodeType } from '@/app/fragments/components/listing'
import { z } from 'zod'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'

const QueryType = z.object({
  spec: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { spec } = QueryType.parse(req.query)
  const data_metanode = krg.getDataNode(spec)
  if (data_metanode) {
    const { spec: id, meta: { icon: _, ...meta }, kind } = data_metanode
    res.status(200).json({
      id,
      kind,
      type: metanodeType(data_metanode),
      meta,
    })
  }
  const process_metanode = krg.getProcessNode(spec)
  if (process_metanode) {
    const { spec: id, meta: { icon: _, ...meta }, kind, inputs, output, story } = process_metanode
    res.status(200).json({
      id,
      kind,
      type: metanodeType(process_metanode),
      meta,
      inputs: dict.init(dict.items(inputs).map(({ key, value }) => ({ key, value: Array.isArray(value) ? [{ id: value.spec }] : { id: value.spec } }))),
      output: { id: output.spec },
      story: story({}),
    })
  }
  throw new NotFoundError()
})
