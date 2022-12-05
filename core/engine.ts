import { Data, Database, Process, Resolved } from '@/core/FPPRG'
import KRG from '@/core/KRG'
import * as dict from '@/utils/dict'

/**
 * This is a minimally viable scg engine -- given a database,
 *  we'll reconcile process outputs when they are created.
 * 
 * Alternatively we could reconcile process outputs whenever they are
 *  requested from the db but aren't already present, this might be more ideal
 *  but might be tricker in a worker served environment.
 * 
 * Thus the way one uses SCG is by adding process nodes to the instance graph
 *  and awaiting the output.
 * 
 */
export default function create_engine(krg: KRG, db: Database) {
  return db.listen(async (table, record) => {
    if (table === 'process') {
      const instanceProcess = record as Process
      try {
        const metaProcess = krg.getProcessNode(instanceProcess.type)
        if (metaProcess === undefined) throw new Error('Unrecognized process')
        const props = {
          data: instanceProcess.data,
          inputs: dict.init(
            dict.items(await instanceProcess.inputs__outputs()).map(({ key, value }) => {
              // propagate errors
              if (value && value.type === 'Error') {
                throw new Error(`Error in ${instanceProcess.type} caused by error in ${metaProcess.inputs[key as string].spec}`)
              }
              return { key, value: value ? value.value : undefined }
          })),
        }
        if ('prompt' in metaProcess) {
          console.log(`Output comes from data`)
          db.upsertResolved(new Resolved(instanceProcess, instanceProcess.data))
        } else {
          const output = metaProcess.output.codec.encode(await metaProcess.resolve(props))
          console.log(`Calling action ${JSON.stringify(metaProcess.spec)} with props ${JSON.stringify(props)} of type ${JSON.stringify(metaProcess.inputs)} to produce ${JSON.stringify(metaProcess.output.spec)}: ${output}`)
          db.upsertResolved(new Resolved(instanceProcess, new Data(metaProcess.output.spec, output)))
        }
      } catch (e) {
        db.upsertResolved(new Resolved(instanceProcess, new Data('Error', JSON.stringify(e))))
      }
    }
  })
}
