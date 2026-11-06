/**
 * Fully Persistent Process Resolution Graph (FPKRG)
 *
 * This consists of:
 * - a fully persistent list data structure for deduplicated lists of processes
 * - process nodes with parent links forming a process graph
 * - resolution nodes attached to each process storing the resolution of that process
 *
 * Ultimately the Process table here correspond to metanode process types,
 *  while the Resolved table corresponds to metanode data types.
 */

import uuid from '@/utils/uuid'
import * as dict from '@/utils/dict'
import * as math from '@/utils/math'
import { z } from 'zod'
import { sql } from 'kysely'
import type { DB, EventBus } from '@/db'
import { NotFoundError, ResponseCodedError, TimeoutError } from '@/spec/error'

/**
 * A process is unique by its own data (configuration) and process outputs
 *  which come before it. In that way, Processes operate like their own
 *  DAG alongside the FPL.
 */
export class Process {
  id: string
  _dependencies: number | undefined = undefined

  constructor(
    /**
     * The type of process, typically corresponding to the KRG Process Type
     */
    public type: string,
    /**
     * Data associated with this process, can used for configuration / prompts
     */
    public data: Data | undefined = undefined,
    /**
     * The input process objects, corresponding to the dependencies of this process
     */
    public inputs: Record<string|number, Process> = {},
    /**
     * The timestamp when this process was saved
     */
    public timestamp: Date | undefined = undefined,
    /**
     * A database instance to operate this class like an ORM
     */
    public db: FPPRG | undefined = undefined,
  ) {
    this.id = uuid([
      type,
      data?.id ?? null,
      dict.sortedItems(inputs).map(({ key, value }) => ({ key, value: value.id })),
      timestamp?.getUTCMilliseconds() ?? null,
    ])
  }

  toJSONWithOutput = async () => {
    const output = await this.output()
    return {
      ...this.toJSON(),
      output: output ? output.toJSON() : null,
    }
  }

  toJSON = () => {
    return {
      'id': this.id,
      'type': this.type,
      'data': this.data?.toJSON() ?? null,
      'timestamp': this.timestamp ? this.timestamp.toISOString() : null,
      'inputs': dict.init(
        dict.sortedItems(this.inputs)
          .map(({ key, value }) => ({ key, value: { id: value.id } }))
      ),
    }
  }

  /**
   * The number of inputs this process depends on,
   *  this is a memoized recursive method across processes.
   * Each process counts itself + the value for all of its inputs.
   */
  dependencies = async (): Promise<number> => {
    if (this._dependencies === undefined) {
      this._dependencies = math.sum(await Promise.all(dict.items(this.inputs).map(({ value }) => value.dependencies()))) + 1
    }
    return this._dependencies
  }

  /**
   * The outputs of all input processes
   */
  inputs__outputs = async (): Promise<Record<string | number | symbol, Data | undefined>> => {
    return dict.init(
      await Promise.all(
        dict.sortedItems(this.inputs).map(async ({ key, value }) =>
          ({ key, value: await value.output() })
        )
      )
    )
  }
  /**
   * The output of this process
   */
  output = async () => {
    if (this.db === undefined) throw new Error('Process not attached to a db')
    return (await this.db.awaitResolved(this.id)).data
  }
}

/**
 * Data resolved by a process
 */
export class Resolved {
  public id: string

  constructor(public process: Process, public data: Data | undefined, public created: Date | undefined = undefined) {
    this.id = process.id
  }

  toJSON = () => {
    return {
      'id': this.id,
      'data': this.data?.toJSON() ?? null,
      'created': this.created ? this.created.toISOString() : null,
    }
  }
}

/**
 * Fully Persistent List of processes to maintain ordering
 */
export class FPL {
  id: string
  constructor(
    public process: Process,
    public parent: FPL | undefined = undefined,
    public cell_metadata: CellMetadata | undefined = undefined,
    public playbook_metadata: PlaybookMetadata | undefined = undefined,
  ) {
    this.id = uuid([
      process.id,
      cell_metadata ?? null,
      playbook_metadata ?? null,
      parent?.id ?? null,
    ])
  }

  toJSONWithOutput = async () => {
    return {
      id: this.id,
      cell_metadata: this.cell_metadata ? this.cell_metadata.toJSON() : null,
      playbook_metadata: this.playbook_metadata ? this.playbook_metadata.toJSON() : null,
      process: await this.process.toJSONWithOutput(),
    }
  }

  toJSON = () => {
    return {
      id: this.id,
      cell_metadata: this.cell_metadata ? this.cell_metadata.toJSON() : null,
      playbook_metadata: this.playbook_metadata ? this.playbook_metadata.toJSON() : null,
      process: this.process.toJSON(),
    }
  }

  /**
   * Resolve the complete list from the FPL data structure
   */
  resolve = () => {
    const fpl: FPL[] = [this]
    let head: FPL = this
    while (head.parent !== undefined) {
      head = head.parent
      fpl.push(head)
    }
    fpl.reverse()
    return fpl
  }

  /**
   * Opposite of resolve, build a fully persistent list out
   *  of a list of processes.
   */
  static fromProcessArray = (processArray: Array<Process>, playbook_metadata = undefined as PlaybookMetadata | undefined) => {
    if (processArray.length === 0) {
      return undefined
    }
    const P = [...processArray]
    const last = P.pop()
    P.reverse()
    let head = undefined
    while (P.length !== 0) {
      head = new FPL(P.pop() as Process, head)
    }
    head = new FPL(last as Process, head, undefined, playbook_metadata)
    return head
  }

  /**
   * Extend an LPL with a new process
   */
  extend = (process: Process) => {
    return new FPL(process, this)
  }

  /**
   * Rebase an LPL's metadata (preserving processes)
   */
  rebasePlaybookMetadata = (old_fpl: FPL, playbook_metadata: PlaybookMetadata) => {
    const fpl: FPL[] = [this]
    let head: FPL | undefined = this
    while (head.parent !== undefined) {
      if (head.id === old_fpl.id) break
      head = head.parent
      fpl.push(head)
    }
    // if we didn't end up on the old_proces, it's not in the FPL
    if (head.id !== old_fpl.id) {
      throw new Error(`${old_fpl.id} not found in FPL`)
    }
    // Walk forward updating processes, replacing any dependencies with the updated one
    head = new FPL(head.process, head.parent, head.cell_metadata, playbook_metadata)
    const rebased = head
    fpl.pop()
    fpl.reverse()
    for (const el of fpl) {
      head = new FPL(el.process, head, el.cell_metadata, el.playbook_metadata)
    }
    return { rebased, head }
  }

  /**
   * Rebase an LPL's cell metadata (preserving processes)
   */
  rebaseCellMetadata = (updates: Record<string, CellMetadata>) => {
    const fplArray = this.resolve()
    let head: FPL | undefined = undefined
    for (const cur of fplArray) {
      if (cur.id in updates) {
        head = new FPL(cur.process, head, updates[cur.id], cur.playbook_metadata)
      } else {
        head = new FPL(cur.process, head, cur.cell_metadata, cur.playbook_metadata)
      }
    }
    return head
  }

  /**
   * Rebase an LPL, this must deal with invalidated processes
   */
  rebase = (old_process: Process, new_process: Process) => {
    // Same as resolve but only up to old_process
    const fpl: FPL[] = [this]
    let head: FPL | undefined = this
    while (head.parent !== undefined) {
      if (head.process.id === old_process.id) break
      head = head.parent
      fpl.push(head)
    }
    // if we didn't end up on the old_proces, it's not in the FPL
    if (head.process.id !== old_process.id) {
      throw new Error(`${old_process} not found in FPL`)
    }
    // Replace old_process with new process
    const update = { [head.process.id]: new_process }
    // Walk forward updating processes, replacing any dependencies with the updated one
    head = new FPL(new_process, head.parent, head.cell_metadata, head.playbook_metadata)
    const rebased = head
    fpl.pop()
    fpl.reverse()
    for (const el of fpl) {
      const new_proc = new Process(
        (update[el.process.id] || el.process).type,
        (update[el.process.id] || el.process).data,
        dict.init(
          dict.sortedItems(el.process.inputs).map(({ key, value }) =>
            ({ key, value: update[value.id] || value })
          )
        ),
        (update[el.process.id] || el.process).timestamp,
        el.process.db,
      )
      head = new FPL(new_proc, head, el.cell_metadata, el.playbook_metadata)
      update[el.process.id] = new_proc
    }
    return { rebased, head }
  }

  rebaseInTime = async () => {
    const fplArray = this.resolve()
    const update: Record<string, Process> = {}
    const now = new Date()
    let head: FPL | undefined = undefined
    for (const el of fplArray) {
      const existing_resolved = el.process.db ? await el.process.db.getResolved((update[el.process.id] || el.process).id) : undefined
      const timestamp = existing_resolved?.created ?? now
      let new_proc = new Process(
        (update[el.process.id] || el.process).type,
        (update[el.process.id] || el.process).data,
        dict.init(
          dict.sortedItems(el.process.inputs).map(({ key, value }) =>
            ({ key, value: update[value.id] || value })
          )
        ),
        timestamp,
        el.process.db,
      )
      if (el.process.db && existing_resolved) {
        new_proc = await el.process.db.upsertProcess(new_proc)
        if (existing_resolved.data) await el.process.db.upsertResolved(new Resolved(new_proc, existing_resolved.data, timestamp))
        await el.process.db.deleteResolved(existing_resolved)
      }
      head = new FPL(new_proc, head, el.cell_metadata, el.playbook_metadata)
      update[el.process.id] = new_proc
    }
    return head
  }

}

/**
 * Content-addressable storage
 */
export class Data {
  id: string

  constructor(public type: string, public value: any) {
    this.id = uuid([type, value])
  }

  toJSON = () => {
    return {
      'id': this.id,
      'type': this.type,
      'value': this.value,
    }
  }
}

/**
 * De-coupled metadata for an FPL
 */
export class CellMetadata {
  id: string

  constructor(
    public label?: string,
    public description?: string,
    public process_visible = true,
    public data_visible = true,
  ) {
    this.id = uuid([
      label,
      description,
      process_visible,
      data_visible,
    ])
  }

  toJSON = () => {
    return {
      'id': this.id,
      'label': this.label,
      'description': this.description,
      'process_visible': this.process_visible,
      'data_visible': this.data_visible,
    }
  }
}

/**
 * De-coupled metadata for an FPL
 */
const SummaryC = z.enum(['auto', 'gpt', 'manual'])

export class PlaybookMetadata {
  id: string

  constructor(
    public title = '',
    public description = '',
    public summary = 'auto' as 'auto' | 'gpt' | 'manual',
    public gpt_summary = '',
  ) {
    this.id = uuid([
      title,
      description,
      summary,
      gpt_summary,
    ])
  }

  toJSON = () => {
    return {
      'id': this.id,
      'title': this.title,
      'description': this.description,
      'summary': this.summary,
      'gpt_summary': this.gpt_summary,
    }
  }
}

export const DataC = z.object({ type: z.string(), value: z.any() })
export const IdOrDataC = z.union([
  z.object({ id: z.string() }),
  DataC,
])
export type IdOrData = z.infer<typeof IdOrDataC>

export const IdOrCellMetadataC = z.union([
  z.object({ id: z.string() }),
  z.object({
    label: z.string().optional(),
    description: z.string().optional(),
    process_visible: z.boolean(),
    data_visible: z.boolean(),
  }),
])
export type IdOrCellMetadata = z.infer<typeof IdOrCellMetadataC>

export const IdOrPlaybookMetadataC = z.union([
  z.object({ id: z.string() }),
  z.object({
    title: z.string().optional(),
    description: z.string().optional(),
    summary: z.enum(['auto', 'gpt', 'manual']).optional(),
    gpt_summary: z.string().optional(),
  }),
])
export type IdOrPlaybookMetadata = z.infer<typeof IdOrPlaybookMetadataC>

export const IdOrProcessC = z.union([
  z.object({
    id: z.string(),
    type: z.string(),
    data: IdOrDataC.optional(),
    inputs: z.record(z.string(), z.object({ id: z.string() })).optional(),
    output: z.object({ id: z.string() }).optional(),
  }),
  z.object({
    id: z.string(),
  }),
  z.object({
    type: z.string(),
    data: IdOrDataC.optional(),
    inputs: z.record(z.string(), z.object({ id: z.string() })).optional(),
    output: z.object({ id: z.string() }).optional(),
  }),
])
export type IdOrProcess = z.infer<typeof IdOrProcessC>

export type DatabaseKeyedTables = { process: Process, fpl: FPL, data: Data, resolved: Resolved }
export type DatabaseListenCallback = <T extends keyof DatabaseKeyedTables>(table: T, record: DatabaseKeyedTables[T]) => void

/**
 * The database maintains records for each table type
 */
export default class FPPRG {
  private processTable: Record<string, Process> = {}
  private fplTable: Record<string, FPL> = {}
  private dataTable: Record<string, Data> = {}
  private cellMetadataTable: Record<string, CellMetadata> = {}
  private playbookMetadataTable: Record<string, PlaybookMetadata> = {}

  constructor(public db: DB, public bus: EventBus) {}

  resolveProcess = async (process: IdOrProcess): Promise<Process> => {
    if ('id' in process) {
      return await this.getProcess(process.id) as Process
    } else {
      return await this.upsertProcess(new Process(
        process.type,
        process.data !== undefined ? await this.resolveData(process.data) : undefined,
        process.inputs !== undefined ? dict.init(await Promise.all(dict.sortedItems(process.inputs).map(async ({ key, value }) => ({ key, value: await this.resolveProcess(value) })))) : {}
      ))
    }
  }
  getProcess = async (id: string): Promise<Process | undefined> => {
    if (!(id in this.processTable)) {
      const result = await this.getProcessRecord(id)
      if (result !== undefined) {
        const process = new Process(
          result.type,
          result.data !== null ? await this.getData(result.data) : undefined,
          dict.init(await Promise.all(dict.sortedItems(result.inputs).map(async ({ key, value }) => ({ key, value: (await this.getProcess(value)) as Process })))),
          result.timestamp ? new Date(result.timestamp) : undefined,
          this
        )
        if (id === process.id) this.processTable[process.id] = process
        else this.processTable[id] = await this.upsertProcess(process)
      }
    }
    return this.processTable[id]
  }
  upsertProcess = async (process: Process) => {
    if (!(process.id in this.processTable)) {
      process.db = this
      // TODO: get these things done in the context manager
      //       in a single transaction
      if (process.data !== undefined) {
        process.data = await this.upsertData(process.data)
      }
      for (const key in process.inputs) {
        const value = process.inputs[key]
        process.inputs[key] = await this.upsertProcess(value)
      }
      await this.db.transaction().execute(async (trx) => {
        await trx.insertInto('process')
          .values({
            id: process.id,
            type: process.type,
            data: process.data !== undefined ? process.data.id : null,
            timestamp: process.timestamp !== undefined ? process.timestamp.toISOString() : null,
          })
          .onConflict(oc => oc.column('id').doUpdateSet({ hits: sql`${sql.ref('process.hits')} + 1` }))
          .execute()
        const inputs = dict.sortedItems(process.inputs).map(({ key, value }) => ({ id: process.id, key: String(key), value: value.id }))
        if (inputs.length > 0) {
          await trx.insertInto('process_input')
            .values(inputs)
            .onConflict(oc => oc.columns(['id', 'key']).doNothing())
            .execute()
        }
      })
      this.processTable[process.id] = process
    } else {
      await this.hit('process', process.id)
    }
    return this.processTable[process.id] as Process
  }

  resolveCellMetadata = async (metadata: IdOrCellMetadata) => {
    if ('id' in metadata) {
      return await this.getCellMetadata(metadata.id) as CellMetadata
    } else {
      return await this.upsertCellMetadata(new CellMetadata(
        metadata.label,
        metadata.description,
        metadata.process_visible,
        metadata.data_visible,
      ))
    }
  }
  getCellMetadata = async (id: string) => {
    if (!(id in this.cellMetadataTable)) {
      const result = await this.db.selectFrom('cell_metadata').selectAll().where('id', '=', id).executeTakeFirst()
      if (result !== undefined) {
        const cellMetadata = new CellMetadata(
          result.label ?? undefined,
          result.description ?? undefined,
          result.process_visible ?? undefined,
          result.data_visible ?? undefined,
        )
        if (id === cellMetadata.id) this.cellMetadataTable[cellMetadata.id] = cellMetadata
        else this.cellMetadataTable[id] = await this.upsertCellMetadata(cellMetadata)
      }
    }
    return this.cellMetadataTable[id] as CellMetadata | undefined
  }
  upsertCellMetadata = async (cell_metadata: CellMetadata) => {
    if (!(cell_metadata.id in this.cellMetadataTable)) {
      await this.db.insertInto('cell_metadata')
        .values({
          id: cell_metadata.id,
          label: cell_metadata.label,
          description: cell_metadata.description,
          process_visible: cell_metadata.process_visible,
          data_visible: cell_metadata.data_visible,
        })
        .onConflict(oc => oc.column('id').doNothing())
        .execute()
      this.cellMetadataTable[cell_metadata.id] = cell_metadata
    }
    return this.cellMetadataTable[cell_metadata.id] as CellMetadata
  }

  resolvePlaybookMetadata = async (playbook_metadata: IdOrPlaybookMetadata) => {
    if ('id' in playbook_metadata) {
      let metadata = await this.getPlaybookMetadata(playbook_metadata.id)
      if (metadata) return metadata
    }
    return await this.upsertPlaybookMetadata(
      new PlaybookMetadata(
        ('title' in playbook_metadata ? playbook_metadata.title : undefined),
        ('description' in playbook_metadata ? playbook_metadata.description : undefined),
        ('summary' in playbook_metadata ? playbook_metadata.summary : undefined) || 'auto',
        ('gpt_summary' in playbook_metadata ? playbook_metadata.gpt_summary : undefined)|| '',
      )
    )
  }
  getPlaybookMetadata = async (id: string) => {
    if (!(id in this.playbookMetadataTable)) {
      const result = await this.db.selectFrom('playbook_metadata').selectAll().where('id', '=', id).executeTakeFirst()
      if (result !== undefined) {
        const playbookMetadata = new PlaybookMetadata(
          result.title ?? undefined,
          result.description ?? undefined,
          SummaryC.optional().catch(undefined).parse(result.summary ?? undefined),
          result.gpt_summary ?? undefined,
        )
        if (id === playbookMetadata.id) this.playbookMetadataTable[playbookMetadata.id] = playbookMetadata
        else this.playbookMetadataTable[id] = await this.upsertPlaybookMetadata(playbookMetadata)
      }
    }
    return this.playbookMetadataTable[id] as PlaybookMetadata | undefined
  }
  upsertPlaybookMetadata = async (playbook_metadata: PlaybookMetadata) => {
    if (!(playbook_metadata.id in this.playbookMetadataTable)) {
      await this.db.insertInto('playbook_metadata')
        .values({
          id: playbook_metadata.id,
          title: playbook_metadata.title,
          description: playbook_metadata.description,
          summary: playbook_metadata.summary,
          gpt_summary: playbook_metadata.gpt_summary,
        })
        .onConflict(oc => oc.column('id').doNothing())
        .execute()
      this.playbookMetadataTable[playbook_metadata.id] = playbook_metadata
    }
    return this.playbookMetadataTable[playbook_metadata.id] as PlaybookMetadata
  }

  resolveData = async (data: IdOrData): Promise<Data> => {
    let resolved: Data | undefined
    if ('id' in data) {
      resolved = await this.getData(data.id)
    }
    if (resolved === undefined && 'type' in data && 'value' in data) {
      resolved = await this.upsertData(new Data(data.type, data.value))
    }
    return resolved as Data
  }
  getData = async (id: string): Promise<Data | undefined> => {
    if (!(id in this.dataTable)) {
      const result = await this.db.selectFrom('data').select(['type', 'value']).where('id', '=', id).executeTakeFirst()
      if (result !== undefined) {
        const data = new Data(
          result.type,
          JSON.parse(result.value),
        )
        if (id === data.id) this.dataTable[data.id] = data
        else this.dataTable[id] = await this.upsertData(data)
      }
    }
    return this.dataTable[id]
  }
  upsertData = async (data: Data): Promise<Data> => {
    if (!(data.id in this.dataTable)) {
      await this.db.insertInto('data')
        .values({ id: data.id, type: data.type, value: JSON.stringify(data.value) })
        .onConflict(oc => oc.column('id').doUpdateSet({ hits: sql`${sql.ref('data.hits')} + 1` }))
        .execute()
      this.dataTable[data.id] = data
    } else {
      await this.hit('data', data.id)
    }
    return this.dataTable[data.id]
  }

  /**
   * The resolution of a process, if it's resolved.
   *  Resolutions can be deleted (to recompute them) so unlike the content-addressed
   *  records, they're always read from the database rather than cached.
   */
  getResolved = async (id: string): Promise<Resolved | undefined> => {
    const result = await this.db.selectFrom('resolved').selectAll().where('id', '=', id).executeTakeFirst()
    if (result === undefined) return undefined
    const process = await this.getProcess(result.id)
    if (process === undefined) return undefined
    const resolved = new Resolved(
      process,
      result.data !== null ? await this.getData(result.data) : undefined,
      result.created !== null ? new Date(result.created) : undefined,
    )
    // the process id changed (e.g. its hashing changed), record it under the new id
    if (id !== resolved.id) await this.upsertResolved(resolved)
    return resolved
  }
  /**
   * Like getResolved but if it's not in the db yet, request & wait for it
   */
  awaitResolved = async (id: string): Promise<Resolved> => {
    const existing = await this.getResolved(id)
    if (existing) return existing
    const process = await this.getProcess(id)
    if (!process) throw new Error('Not found')
    await new Promise<void>(async (resolve, reject) => {
      const ctx = { resolved: false }
      const unsub = this.bus.listen(`insert:resolved`, (data) => {
        const record = z.object({ id: z.string() }).parse(data)
        if (record.id === id) {
          if (!ctx.resolved) {
            ctx.resolved = true
            unsub()
            resolve()
          }
        }
      })
      console.debug(`awaiting ${process.type} (${id})`)
      await this.bus.send('work-queue', { id, priority: - (await process.dependencies()) })
      setTimeout(() => {
        if (!ctx.resolved) {
          ctx.resolved = true
          unsub()
          reject(new TimeoutError())
        }
      }, 5000)
    })
    const resolved = await this.getResolved(id)
    if (resolved === undefined) throw new Error(`Received event but not resolved (${id})`)
    return resolved
  }
  upsertResolved = async (resolved: Resolved): Promise<Resolved> => {
    if (resolved.data) {
      await this.upsertData(resolved.data)
    }
    const inserted = await this.db.insertInto('resolved')
      .values({
        id: resolved.id,
        data: resolved.data !== undefined ? resolved.data.id : null,
        created: (resolved.created ?? new Date()).toISOString(),
      })
      .onConflict(oc => oc.column('id').doNothing())
      .executeTakeFirst()
    if (Number(inserted.numInsertedOrUpdatedRows ?? 0) > 0) this.bus.inserted('resolved', resolved.id)
    return resolved
  }
  deleteResolved = async (resolved: Resolved) => {
    await this.db.deleteFrom('resolved').where('id', '=', resolved.id).execute()
  }

  resolveFPL = async ({ data = {}, workflow, metadata }: { data?: Record<string, z.infer<typeof DataC>>, workflow: z.infer<typeof IdOrProcessC>[], metadata?: z.infer<typeof IdOrPlaybookMetadataC> }) => {
    // resolve the process array by walking through the specification, collecting instantiated processes
    const processArray: Process[] = []
    const processArrayLookup: Record<string|number, string> = {}
    for (const k in data) {
      await this.resolveData(data[k])
    }
    for (const el of workflow) {
      if ('id' in el && 'type' in el) {
        // given both id and type makes this a 
        const { id, ...proc } = el
        if (proc.data !== undefined) {
          if ('id' in proc.data) {
            proc.data = data[proc.data.id]
          }
        }
        if (proc.inputs !== undefined) {
          for (const k in proc.inputs) {
            if (!(proc.inputs[k].id in processArrayLookup)) throw new ResponseCodedError(400, `${proc.inputs[k].id} not found in preceeding graph`)
            proc.inputs[k].id = processArrayLookup[proc.inputs[k].id]
          }
        }
        const resolvedProc = await this.resolveProcess(proc)
        processArrayLookup[id] = resolvedProc.id
        processArray.push(resolvedProc)
        if (proc.output !== undefined) {
          await this.upsertResolved(new Resolved(resolvedProc, await this.resolveData(proc.output), resolvedProc.timestamp ?? new Date()))
        }
      } else if ('id' in el) {
        const resolvedProc = await this.resolveProcess(el)
        processArrayLookup[resolvedProc.id] = resolvedProc.id
        processArray.push(resolvedProc)
      } else {
        const resolvedProc = await this.resolveProcess(el)
        processArrayLookup[resolvedProc.id] = resolvedProc.id
        processArray.push(resolvedProc)
      }
    }
    const playbookMetadata = metadata ? await this.resolvePlaybookMetadata(metadata) : undefined
    const processArrayFPL = FPL.fromProcessArray(processArray, playbookMetadata)
    if (!processArrayFPL) throw new NotFoundError()
    return await this.upsertFPL(processArrayFPL)
  }
  getFPL = async (id: string): Promise<FPL | undefined> => {
    if (!(id in this.fplTable)) {
      const result = await this.db.selectFrom('fpl').selectAll().where('id', '=', id).executeTakeFirst()
      if (result !== undefined) {
        const fpl = new FPL(
          (await this.getProcess(result.process)) as Process,
          result.parent !== null ? (await this.getFPL(result.parent)) as FPL : undefined,
          result.cell_metadata ? await this.getCellMetadata(result.cell_metadata) : undefined,
          result.playbook_metadata ? await this.getPlaybookMetadata(result.playbook_metadata) : undefined,
        )
        if (id === fpl.id) this.fplTable[fpl.id] = fpl
        else this.fplTable[id] = await this.upsertFPL(fpl)
      }
    }
    return this.fplTable[id]
  }
  upsertFPL = async (fpl: FPL): Promise<FPL> => {
    if (!(fpl.id in this.fplTable)) {
      fpl.process = await this.upsertProcess(fpl.process)
      if (fpl.parent !== undefined) {
        fpl.parent = await this.upsertFPL(fpl.parent)
      }
      if (fpl.cell_metadata !== undefined) {
        fpl.cell_metadata = await this.upsertCellMetadata(fpl.cell_metadata)
      }
      if (fpl.playbook_metadata !== undefined) {
        fpl.playbook_metadata = await this.upsertPlaybookMetadata(fpl.playbook_metadata)
      }
      try{
        await this.db.insertInto('fpl')
          .values({
            id: fpl.id,
            process: fpl.process.id,
            cell_metadata: fpl.cell_metadata !== undefined ? fpl.cell_metadata.id : null,
            playbook_metadata: fpl.playbook_metadata !== undefined ? fpl.playbook_metadata.id : null,
            parent: fpl.parent !== undefined ? fpl.parent.id : null,
          })
          .onConflict(oc => oc.column('id').doUpdateSet({ hits: sql`${sql.ref('fpl.hits')} + 1` }))
          .execute()
      }catch(e){
        console.error({ fpl_id: fpl.id, fpl_process_id: fpl.process.id, fpl_parent_id: fpl.parent?.id })
        throw e
      }
      this.fplTable[fpl.id] = fpl
    } else {
      await this.hit('fpl', fpl.id)
    }
    return this.fplTable[fpl.id]
  }

  /**
   * A process along with its inputs
   */
  private getProcessRecord = async (id: string) => {
    const process = await this.db.selectFrom('process')
      .select(['id', 'type', 'data', 'timestamp'])
      .where('id', '=', id)
      .executeTakeFirst()
    if (process === undefined) return undefined
    const inputs = await this.db.selectFrom('process_input')
      .select(['key', 'value'])
      .where('id', '=', id)
      .execute()
    return {
      id: process.id,
      type: process.type,
      data: process.data,
      timestamp: process.timestamp,
      inputs: dict.init(inputs.map(({ key, value }) => ({ key, value }))),
    }
  }

  /**
   * Count a usage of an existing record
   */
  private hit = async (table: 'process' | 'data' | 'fpl', id: string) => {
    await this.db.updateTable(table)
      .set({ hits: sql`hits + 1` })
      .where('id', '=', id)
      .execute()
  }

  dump = async () => {
    return {
      dataTable: this.dataTable,
      processTable: this.processTable,
      fplTable: this.fplTable,
      cellMetadataTable: this.playbookMetadataTable,
      playbookMetadataTable: this.cellMetadataTable,
    }
  }
}
