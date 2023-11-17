import fpprg from '@/app/fpprg'
import { z } from 'zod'
import handler from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import * as dict from '@/utils/dict'

const QueryType = z.object({
  fpl_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError()
  const { fpl_id } = QueryType.parse(req.query)
  const fpl = await fpprg.getFPL(fpl_id)
  if (fpl === undefined) throw new NotFoundError()
  const data: Record<string, any> = {}
  const workflow = []
  for (const el of fpl.resolve()) {
    workflow.push({
      id: el.process.id,
      type: el.process.type,
      inputs: !dict.isEmpty(el.process.inputs) ? dict.init(dict.items(el.process.inputs).map(({ key, value }) => ({ key, value: { id: value.id } }))) : undefined,
      data: el.process.data ? { id: el.process.data.id } : undefined,
    })
    if (el.process.data) {
      const { id: _, ...data_node } = el.process.data.toJSON()
      data[el.process.data.id] = data_node
    }
  }
  res.status(200).json({
    data,
    workflow,
    metadata: fpl.playbook_metadata ? fpl.playbook_metadata.toJSON() : undefined
  })
})
