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
import { Db } from '@/utils/orm'
import * as fpprgSchema from '@/db/fpprg'
import { TypedSchema } from '@/spec/sql'
import { TimeoutError } from '@/spec/error'

/**
 * A process is unique by its own data (configuration) and process outputs
 *  which come before it. In that way, Processes operate like their own
 *  DAG alongside the FPL.
 */
export class Process {
  id: string
  resolved: Resolved | undefined = undefined
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
     * A database instance to operate this class like an ORM
     */
    public db: FPPRG | undefined = undefined,
  ) {
    this.id = uuid([type, data?.id, dict.sortedItems(inputs).map(({ key, value }) => ({ key, value: value.id }))])
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
      'data': this.data !== undefined ? this.data.toJSON() : null,
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
    if (this.resolved !== undefined) return this.resolved.data
    if (this.db === undefined) throw new Error('Process not attached to a db')
    this.resolved = await this.db.awaitResolved(this.id)
    return this.resolved.data
  }
}

/**
 * Data resolved by a process
 */
export class Resolved {
  public id: string

  constructor(public process: Process, public data: Data | undefined) {
    this.id = process.id
  }

  toJSON = () => {
    return {
      'id': this.id,
      'data': this.data ? this.data.toJSON() : null,
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
      cell_metadata,
      playbook_metadata,
      parent?.id,
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
    let head = undefined
    for (let i = 0; i < fplArray.length; i++) {
      const cur = fplArray[i]
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
        el.process.db,
      )
      head = new FPL(new_proc, head, el.cell_metadata, el.playbook_metadata)
      update[el.process.id] = new_proc
    }
    return { rebased, head }
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
  }),
  z.object({
    id: z.string(),
  }),
  z.object({
    type: z.string(),
    data: IdOrDataC.optional(),
    inputs: z.record(z.string(), z.object({ id: z.string() })).optional(),
  }),
])
export type IdOrProcess = z.infer<typeof IdOrProcessC>

export type DatabaseKeyedTables = { process: Process, fpl: FPL, data: Data, resolved: Resolved }
export type DatabaseListenCallback = <T extends keyof DatabaseKeyedTables>(table: T, record: DatabaseKeyedTables[T]) => void

type Schema = Pick<typeof fpprgSchema, 'data' | 'process' | 'process_input' | 'resolved' | 'process_complete' | 'fpl' | 'cell_metadata' | 'playbook_metadata'>
type TofSchema<T> = T extends TypedSchema<infer T_> ? T_ : never
type DbSchemaT = { [K in keyof Schema]: TofSchema<Schema[K]> }

/**
 * The database maintains records for each table type
 */
export default class FPPRG {
  private processTable: Record<string, Process> = {}
  private fplTable: Record<string, FPL> = {}
  private resolvedTable: Record<string, Resolved> = {}
  private dataTable: Record<string, Data> = {}
  private cellMetadataTable: Record<string, CellMetadata> = {}
  private playbookMetadataTable: Record<string, PlaybookMetadata> = {}

  constructor(public db: Db<DbSchemaT>) {}

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
  getProcess = async (id: string, force = false): Promise<Process | undefined> => {
    if (!(id in this.processTable) || force) {
      const result = await this.db.objects.process_complete.findUnique({ where: { id } })
      if (result !== null) {
        const process = new Process(
          result.type,
          result.data !== null ? await this.getData(result.data) : undefined,
          dict.init(await Promise.all(dict.sortedItems(result.inputs).map(async ({ key, value }) => ({ key, value: (await this.getProcess(value)) as Process })))),
          this
        )
        if (id === process.id) this.processTable[process.id] = process
        else this.processTable[id] = await this.upsertProcess(process)
        // this process was already resolved?
        // save that information, saves a lookup or two
        if (result.resolved) {
          if (!(id in this.resolvedTable)) {
            const resolved = new Resolved(
              process,
              result.output !== null ? (await this.getData(result.output)) : undefined,
            )
            if (id === resolved.id) this.resolvedTable[resolved.id] = process.resolved = resolved
            else this.resolvedTable[id] = process.resolved = await this.upsertResolved(resolved)
          }
        }
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
      // TODO:
      // await this.db.objects.process.upsert({
      //   where: { id: process.id },
      //   update: {
      //     type: process.type,
      //     data: process.data !== undefined ? process.data.id : null,
      //   },
      //   create: {
      //     id: process.id,
      //     type: process.type,
      //     data: process.data !== undefined ? process.data.id : null,
      //     inputs: dict.sortedItems(process.inputs).map(({ key, value }) => ({
      //       id: process.id,
      //       key,
      //       value: value.id,
      //     }))
      //   }
      // })
      await this.db.objects.process.upsert({
        where: { id: process.id },
        create: {
          id: process.id,
          type: process.type,
          data: process.data !== undefined ? process.data.id : null,
        }
      })
      for (const key in process.inputs) {
        const value = process.inputs[key]
        await this.db.objects.process_input.upsert({
          where: { id: process.id, key },
          create: {
            id: process.id,
            key,
            value: value.id,
          }
        })
      }
      this.processTable[process.id] = process
    } else {
      await this.db.objects.process.update({
        where: {
          id: process.id
        },
        data: {
          id: process.id
        }
      })
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
      const result = await this.db.objects.cell_metadata.findUnique({ where: { id } })
      if (result !== null) {
        const cellMetadata = new CellMetadata(
          result.label,
          result.description,
          result.process_visible,
          result.data_visible,
        )
        if (id === cellMetadata.id) this.cellMetadataTable[cellMetadata.id] = cellMetadata
        else this.cellMetadataTable[id] = await this.upsertCellMetadata(cellMetadata)
      }
    }
    return this.cellMetadataTable[id] as CellMetadata | undefined
  }
  upsertCellMetadata = async (cell_metadata: CellMetadata) => {
    if (!(cell_metadata.id in this.cellMetadataTable)) {
      await this.db.objects.cell_metadata.upsert({
        where: { id: cell_metadata.id },
        create: {
          id: cell_metadata.id,
          label: cell_metadata.label,
          description: cell_metadata.description,
          process_visible: cell_metadata.process_visible,
          data_visible: cell_metadata.data_visible,
        }
      })
      this.cellMetadataTable[cell_metadata.id] = cell_metadata
    }
    return this.cellMetadataTable[cell_metadata.id] as CellMetadata
  }

  resolvePlaybookMetadata = async (playbook_metadata: IdOrPlaybookMetadata) => {
    if ('id' in playbook_metadata) {
      return await this.getPlaybookMetadata(playbook_metadata.id) as PlaybookMetadata
    } else {
      return await this.upsertPlaybookMetadata(
        new PlaybookMetadata(
          playbook_metadata.title,
          playbook_metadata.description,
          playbook_metadata.summary || 'auto',
          playbook_metadata.gpt_summary || '',
        ))
    }
  }
  getPlaybookMetadata = async (id: string) => {
    if (!(id in this.playbookMetadataTable)) {
      const result = await this.db.objects.playbook_metadata.findUnique({ where: { id } })
      if (result !== null) {
        const playbookMetadata = new PlaybookMetadata(
          result.title,
          result.description,
          result.summary,
          result.gpt_summary,
        )
        if (id === playbookMetadata.id) this.playbookMetadataTable[playbookMetadata.id] = playbookMetadata
        else this.playbookMetadataTable[id] = await this.upsertPlaybookMetadata(playbookMetadata)
      }
    }
    return this.playbookMetadataTable[id] as PlaybookMetadata | undefined
  }
  upsertPlaybookMetadata = async (playbook_metadata: PlaybookMetadata) => {
    if (!(playbook_metadata.id in this.playbookMetadataTable)) {
      await this.db.objects.playbook_metadata.upsert({
        where: { id: playbook_metadata.id },
        create: {
          id: playbook_metadata.id,
          title: playbook_metadata.title,
          description: playbook_metadata.description,
          summary: playbook_metadata.summary,
          gpt_summary: playbook_metadata.gpt_summary,
        }
      })
      this.playbookMetadataTable[playbook_metadata.id] = playbook_metadata
    }
    return this.playbookMetadataTable[playbook_metadata.id] as PlaybookMetadata
  }

  resolveData = async (data: IdOrData): Promise<Data> => {
    if ('id' in data) {
      return await this.getData(data.id) as Data
    } else {
      return await this.upsertData(new Data(data.type, data.value))
    }
  }
  getData = async (id: string): Promise<Data | undefined> => {
    if (!(id in this.dataTable)) {
      const result = await this.db.objects.data.findUnique({ where: { id } })
      if (result !== null) {
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
      await this.db.objects.data.upsert({
        where: { id: data.id },
        create: {
          id: data.id,
          type: data.type,
          value: JSON.stringify(data.value),
        }
      })
      this.dataTable[data.id] = data
    } else {
      await this.db.objects.data.update({
        where: {
          id: data.id
        },
        data: {
          id: data.id
        }
      })
    }
    return this.dataTable[data.id]
  }

  getResolved = async (id: string): Promise<Resolved | undefined> => {
    if (!(id in this.resolvedTable)) {
      const result = await this.db.objects.resolved.findUnique({ where: { id } })
      if (result !== null) {
        const process = (await this.getProcess(result.id)) as Process
        const resolved = new Resolved(
          process,
          result.data !== null ? (await this.getData(result.data)) as Data : undefined,
        )
        if (id === resolved.id) this.resolvedTable[resolved.id] = process.resolved = resolved
        else this.resolvedTable[id] = process.resolved = await this.upsertResolved(resolved)
      }
    }
    return this.resolvedTable[id]
  }
  /**
   * Like getResolved but if it's not in the db yet, request & wait for it
   */
  awaitResolved = async (id: string): Promise<Resolved> => {
    const resolved = await this.getResolved(id)
    if (!resolved) {
      const process = await this.getProcess(id)
      if (!process) throw new Error('Not found')
      await new Promise<void>(async (resolve, reject) => {
        const ctx = { resolved: false }
        const unsub = this.db.listen((evt, data) => {
          if (evt === `insert:resolved`) {
            const record = z.object({ id: z.string() }).parse(data)
            if (record.id === id) {
              if (!ctx.resolved) {
                ctx.resolved = true
                unsub()
                resolve()
              }
            }
          }
        })
        console.debug(`Awaiting ${process.type} (${id})`)
        await this.db.send('work-queue', { id, priority: - (await process.dependencies()) })
        setTimeout(() => {
          if (!ctx.resolved) {
            ctx.resolved = true
            unsub()
            reject(new TimeoutError())
          }
        }, 5000)
      })
      if (await this.getResolved(id) === undefined) {
        throw new Error(`Received event but not resolved (${id})`)
      }
    }
    return this.resolvedTable[id]
  }
  upsertResolved = async (resolved: Resolved): Promise<Resolved> => {
    if (!(resolved.id in this.resolvedTable)) {
      if (resolved.data) {
        await this.upsertData(resolved.data)
      }
      await this.db.objects.resolved.upsert({
        where: { id: resolved.id },
        create: {
          id: resolved.id,
          data: resolved.data !== undefined ? resolved.data.id : null,
        },
      })
      this.resolvedTable[resolved.id] = resolved
      resolved.process.resolved = resolved
    }
    return this.resolvedTable[resolved.id]
  }
  deleteResolved = async (resolved: Resolved) => {
    if (resolved.id in this.resolvedTable) {
      await this.db.objects.resolved.delete({ where: { id: resolved.id } })
      if (resolved.id in this.processTable) {
        this.processTable[resolved.id].resolved = undefined
      }
      delete this.resolvedTable[resolved.id]
    }
  }

  getFPL = async (id: string): Promise<FPL | undefined> => {
    if (!(id in this.fplTable)) {
      const result = await this.db.objects.fpl.findUnique({ where: { id } })
      if (result !== null) {
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
        await this.db.objects.fpl.upsert({
          where: { id: fpl.id },
          create: {
            id: fpl.id,
            process: fpl.process.id,
            cell_metadata: fpl.cell_metadata !== undefined ? fpl.cell_metadata.id : null,
            playbook_metadata: fpl.playbook_metadata !== undefined ? fpl.playbook_metadata.id : null,
            parent: fpl.parent !== undefined ? fpl.parent.id : null,
          }
        })
      }catch(e){
        console.error({ fpl_id: fpl.id, fpl_process_id: fpl.process.id, fpl_parent_id: fpl.parent?.id })
        throw e
      }
      this.fplTable[fpl.id] = fpl
    } else {
      await this.db.objects.fpl.update({
        where: {
          id: fpl.id
        },
        data: {
          id: fpl.id
        }
      })
    }
    return this.fplTable[fpl.id]
  }

  dump = async () => {
    return {
      dataTable: this.dataTable,
      processTable: this.processTable,
      resolvedTable: this.resolvedTable,
      fplTable: this.fplTable,
      cellMetadataTable: this.playbookMetadataTable,
      playbookMetadataTable: this.cellMetadataTable,
    }
  }
}
