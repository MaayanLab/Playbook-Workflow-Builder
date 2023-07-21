import fpprg from '@/app/fpprg'
import { FPL, Process } from '@/core/FPPRG'
import { z } from 'zod'
import { IdOrProcessC } from '@/core/FPPRG'
import { NotFoundError, ResponseCodedError, UnsupportedMethodError } from '@/spec/error'
import handler from '@/utils/next-rest'

const BodyType = z.array(IdOrProcessC)

export default handler(async (req, res) => {
  if (req.method !== 'POST') throw new UnsupportedMethodError()
  // load the submitted processArray specification
  const processArraySpec = BodyType.parse(req.body)
  // resolve the process array by walking through the specification, collecting instantiated processes
  const processArray: Process[] = []
  const processArrayLookup: Record<string|number, string> = {}
  for (const el of processArraySpec) {
    if ('id' in el && 'type' in el) {
      // given both id and type makes this a 
      const { id, ...proc } = el
      if (proc.inputs !== undefined) {
        for (const k in proc.inputs) {
          if (!(proc.inputs[k].id in processArrayLookup)) throw new ResponseCodedError(400, `${proc.inputs[k].id} not found in preceeding graph`)
          proc.inputs[k].id = processArrayLookup[proc.inputs[k].id]
        }
      }
      const resolvedProc = await fpprg.resolveProcess(proc)
      processArrayLookup[id] = resolvedProc.id
      processArray.push(resolvedProc)
    } else if ('id' in el) {
      const resolvedProc = await fpprg.resolveProcess(el)
      processArrayLookup[resolvedProc.id] = resolvedProc.id
      processArray.push(resolvedProc)
    } else {
      const resolvedProc = await fpprg.resolveProcess(el)
      processArrayLookup[resolvedProc.id] = resolvedProc.id
      processArray.push(resolvedProc)
    }
  }
  const processArrayFPL = FPL.fromProcessArray(processArray)
  if (!processArrayFPL) throw new NotFoundError()
  const fpl = await fpprg.upsertFPL(processArrayFPL)
  res.status(200).json(fpl.id)
})
