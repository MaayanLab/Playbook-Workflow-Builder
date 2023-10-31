import FPPRG, { Data, Process, Resolved } from '@/core/FPPRG'
import KRG from '@/core/KRG'
import * as dict from '@/utils/dict'
import { z } from 'zod'
import * as array from '@/utils/array'
import { UnboundError, TimeoutError } from '@/spec/error'
import { StatusUpdate } from '@/spec/metanode'

const JobC = z.object({ data: z.object({ id: z.string() }) })

export async function decode_complete_process_inputs(krg: KRG, instanceProcess: Process) {
  const metaProcess = krg.getProcessNode(instanceProcess.type)
  if (metaProcess === undefined) throw new Error('Unrecognized process')
  return dict.items(await instanceProcess.inputs__outputs()).reduce<Record<string, unknown>>((inputs, { key, value }) => {
    const [arg, ...i] = (key as string).split(':')
    if (i.length > 0 && !Array.isArray(metaProcess.inputs[arg])) {
      throw new Error('Received multiple args, but not an array')
    } else if (i.length === 0 && Array.isArray(metaProcess.inputs[arg])) {
      throw new Error('Expected multiple args')
    }
    const metaProcessInput = array.ensureOne(metaProcess.inputs[arg])
    if (value === undefined) {
      // handle nodes
      throw new UnboundError()
    } else if (value.type === 'Error') {
      // propagate errors
      throw new Error(`${instanceProcess.type} can't run because of error in ${metaProcessInput.spec}`)
    }
    const value_decoded = value ? metaProcessInput.codec.decode(value.value) : undefined
    if (i.length > 0) {
      return {...inputs, [arg]: [...(inputs[arg]||[]) as unknown[], value_decoded] }
    } else {
      return {...inputs, [arg]: value_decoded }
    }
  }, {})
}
export async function decode_complete_process_output(krg: KRG, instanceProcess: Process) {
  const output = await instanceProcess.output()
  if (!output) throw new Error('No output')
  const metaDataNode = krg.getDataNode(output.type)
  return metaDataNode.codec.decode(output.value)
}

/**
 * Given an instanceProcess, "resolve" its output, this is typically done
 *  by executing it, though prompts are implicitly identity functions.
 */
export async function resolve_process(krg: KRG, instanceProcess: Process, notify: (update: StatusUpdate) => void) {
  try {
    const metaProcess = krg.getProcessNode(instanceProcess.type)
    if (metaProcess === undefined) throw new Error('Unrecognized process')
    console.debug(`Preparing ${metaProcess.spec}`)
    const props = {
      data: instanceProcess.data?.value,
      inputs: await decode_complete_process_inputs(krg, instanceProcess),
      notify,
    }
    console.debug(`Processing ${metaProcess.spec}`)
    if ('prompt' in metaProcess && props.data === undefined) throw new UnboundError()
    const output = metaProcess.output.codec.encode(await metaProcess.resolve(props))
    return new Resolved(instanceProcess, new Data(metaProcess.output.spec, output))
  } catch (e) {
    if (e instanceof TimeoutError) {
      console.warn(e)
      throw e
    } else if (e instanceof UnboundError) {
      return new Resolved(instanceProcess, undefined)
    } else {
      console.error(e)
      return new Resolved(instanceProcess, new Data('Error', (e as Error).toString()))
    }
  }
}

/**
 * This worker receives jobs from the boss work-queue which ensures
 *  submitted jobs are given to one and only one worker. This should
 *  be run with several replicas in production.
 */
export function start_workers(krg: KRG, db: FPPRG, n_workers: number) {
  return db.db.work('work-queue', { teamSize: n_workers, teamConcurrency: n_workers }, async (job) => {
    // the job.data should contain the process id
    const { data: { id: processId } } = JobC.parse(job)
    console.debug(`checking ${processId}..`)
    // we fetch it from the db
    const instanceProcess = await db.getProcess(processId, true)
    if (!instanceProcess) {
      console.error(`process '${processId}' not found`)
    } else if (instanceProcess.resolved === undefined) {
      console.debug(`resolving ${instanceProcess.type} (${processId})..`)
      const resolved = await resolve_process(krg, instanceProcess, (update) => {
        // broadcast partial updates through the db
        db.db.notify(`distributed:resolved:${instanceProcess.id}:status`, JSON.stringify(update))
      })
      // store the result in the db
      await db.upsertResolved(resolved)
      console.debug(`completed ${instanceProcess.type} ${processId}`)
    }
  })
}
