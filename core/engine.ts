import FPPRG, { Data, Process, Resolved, TimeoutError } from '@/core/FPPRG'
import KRG from '@/core/KRG'
import * as dict from '@/utils/dict'
import { z } from 'zod'

const JobC = z.object({ data: z.object({ id: z.string() }) })

/**
 * This error occurs when the input node is not populated yet
 */
class UnboundError extends Error {}

/**
 * Given an instanceProcess, "resolve" its output, this is typically done
 *  by executing it, though prompts are implicitly identity functions.
 */
export async function resolve_process(krg: KRG, instanceProcess: Process) {
  try {
    const metaProcess = krg.getProcessNode(instanceProcess.type)
    if (metaProcess === undefined) throw new Error('Unrecognized process')
    const props = {
      data: instanceProcess.data,
      inputs: dict.init(
        dict.items(await instanceProcess.inputs__outputs()).map(({ key, value }) => {
          if (value === undefined) {
            // handle nodes
            throw new UnboundError()
          } else if (value.type === 'Error') {
            // propagate errors
            throw new Error(`${instanceProcess.type} can't run because of error in ${metaProcess.inputs[key as string].spec}`)
          }
          return { key, value: value ? metaProcess.inputs[key as string].codec.decode(value.value) : undefined }
      })),
    }
    if ('prompt' in metaProcess) {
      console.debug(`Output comes from data`)
      return new Resolved(instanceProcess, instanceProcess.data)
    } else {
      // TODO: add a timeout here
      const output = metaProcess.output.codec.encode(await metaProcess.resolve(props))
      console.debug(`Calling action ${JSON.stringify(metaProcess.spec)} with props ${JSON.stringify(props)} of type ${JSON.stringify(metaProcess.inputs)}`)
      return new Resolved(instanceProcess, new Data(metaProcess.output.spec, output))
    }
  } catch (e) {
    if (e instanceof TimeoutError) {
      throw e
    } else if (e instanceof UnboundError) {
      return new Resolved(instanceProcess, undefined)
    } else {
      console.error(e)
      return new Resolved(instanceProcess, new Data('Error', JSON.stringify((e as Error).toString())))
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
    if (!instanceProcess) throw new Error(`Process '${processId}' not found`)
    if (instanceProcess.resolved === undefined) {
      console.debug(`resolving ${processId}..`)
      const resolved = await resolve_process(krg, instanceProcess)
      // store the result in the db
      await db.upsertResolved(resolved)
      console.debug(`completed ${processId}`)
    }
  })
}
