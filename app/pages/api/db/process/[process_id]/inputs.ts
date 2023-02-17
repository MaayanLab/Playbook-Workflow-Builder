import fpprg from '@/app/fpprg'
import { z } from 'zod'
import handler from '@/utils/next-rest'
import { NotFoundError, UnboundError } from '@/spec/error'
import * as dict from '@/utils/dict'

const QueryType = z.object({
  process_id: z.string(),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new Error('Unsupported method')
  const { process_id } = QueryType.parse(req.query)
  const process = await fpprg.getProcess(process_id)
  if (process === undefined) throw new NotFoundError()
  const inputs = dict.items(await process.inputs__outputs()).reduce<Record<string, any>>((inputs, { key, value }) => {
    const [arg, ...i] = (key as string).split(':')
    if (value === undefined) {
      // handle nodes
      throw new UnboundError()
    } else if (value.type === 'Error') {
      // propagate errors
      throw new Error(`${process.type} can't run because of error in prerequisite`)
    }
    const value_decoded = value ? value.value : undefined
    if (i.length > 0) {
      return {...inputs, [arg]: [...(inputs[arg]||[]), value_decoded] }
    } else {
      return {...inputs, [arg]: value_decoded }
      }
  }, {})
  res.status(200).json(inputs)
})
