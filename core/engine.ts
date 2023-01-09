import { Data, Database, Process, Resolved } from '@/core/FPPRG'
import KRG from '@/core/KRG'
import * as dict from '@/utils/dict'
import { PgDatabase } from '@/core/FPPRG'
import * as z from 'zod'

const JobC = z.object({ id: z.string() })

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
      const output = metaProcess.output.codec.encode(await metaProcess.resolve(props))
      console.debug(`Calling action ${JSON.stringify(metaProcess.spec)} with props ${JSON.stringify(props)} of type ${JSON.stringify(metaProcess.inputs)}`)
      return new Resolved(instanceProcess, new Data(metaProcess.output.spec, output))
    }
  } catch (e) {
    if (e instanceof UnboundError) {
      return new Resolved(instanceProcess, undefined)
    } else {
      console.error(e)
      return new Resolved(instanceProcess, new Data('Error', JSON.stringify((e as Error).toString())))
    }
  }
}

/**
 * This is a minimally viable scg engine -- given a database,
 *  we'll reconcile process outputs when they are created.
 */
export function process_insertion_dispatch(krg: KRG, db: Database) {
  return db.listen(async (table, record) => {
    if (table === 'process') {
      const instanceProcess = record as Process
      const resolved = await resolve_process(krg, instanceProcess)
      await db.upsertResolved(resolved)
    }
  })
}

/**
 * This worker receives jobs from the boss work-queue which ensures
 *  submitted jobs are given to one and only one worker. This should
 *  be run with several replicas in production.
 */
export function start_workers(krg: KRG, db: PgDatabase, n_workers: number) {
  db.boss.on('error', error => console.error(error))
  ;(async () => {
    console.log(`starting worker(s)..`)
    await db.boss.start()
    await db.boss.work('work-queue', { teamSize: n_workers, teamConcurrency: n_workers }, async (job) => {
      // the job.data should contain the process id
      const { id: processId } = JobC.parse(job.data)
      console.debug(`checking ${processId}..`)
      // we fetch it from the db
      const instanceProcess = await db.getProcess(processId, true)
      if (!instanceProcess) throw new Error(`Process ${job.data} not found`)
      if (instanceProcess.resolved === undefined) {
        console.debug(`resolving ${processId}..`)
        const resolved = await resolve_process(krg, instanceProcess)
        // store the result in the db
        await db.upsertResolved(resolved)
        console.debug(`completed ${processId}`)
      }
    })
  })().catch(error => console.error(error))
  return () => {
    db.boss.stop().catch(error => console.error(error))
  }
}
