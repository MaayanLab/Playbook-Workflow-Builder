/**
 * Fully Persistent Process Resolution Graph (FPKRG)
 * 
 * This is an in-memory database implementation of a fully persistent process resolution graph.
 * This consists of:
 * - a fully persistent list data structure for deduplicated lists of processes
 * - process nodes with parent links forming a process graph
 * - resolution nodes attached to each process storing the resolution of that process
 * 
 * Ultimately the Process table here correspond to metanode process types,
 *  while the Resolved table corresponds to metanode data types.
 */

import sha256 from '@/utils/sha256'
import * as dict from '@/utils/dict'
import { z } from 'zod'

/**
 * Content-addressable storage
 */
export class Data {
  id: string

  constructor(public type: string, public value: string) {
    this.id = sha256([type, value])
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
 * A process is unique by its own data (configuration) and process outputs
 *  which come before it. In that way, Processes operate like their own
 *  DAG alongside the FPL.
 */
export class Process {
  id: string

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
    public db: Database | undefined = undefined,
  ) {
    this.id = sha256([type, data, dict.items(inputs).map(({ key, value }) => ({ key, value: value.id }))])
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
        dict.items(this.inputs)
          .map(({ key, value }) => ({ key, value: { id: value.id } }))
      ),
    }
  }

  /**
   * The outputs of all input processes
   */
  inputs__outputs = async (): Promise<Record<string | number | symbol, Data | undefined>> => {
    return dict.init<Data | undefined>(
      await Promise.all(
        dict.items(this.inputs).map(async ({ key, value }) =>
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
    const resolved = await this.db.awaitResolved(this.id)
    return resolved.data
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
  constructor(public process: Process, public parent: FPL | undefined = undefined) {
    this.id = sha256([process.id, parent ? parent.id : null])
  }

  toJSONWithOutput = async () => {
    return {
      id: this.id,
      process: await this.process.toJSONWithOutput(),
    }
  }

  toJSON = () => {
    return {
      'id': this.id,
      'process': this.process.toJSON(),
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
  static fromProcessArray = (processArray: Array<Process>) => {
    const P = [...processArray]
    P.reverse()
    let head = undefined
    while (P.length !== 0) {
      head = new FPL(P.pop() as Process, head)
    }
    return head
  }

  /**
   * Extend an LPL with a new process
   */
  extend = (process: Process) => {
    return new FPL(process, this)
  }
  /**
   * Rebase an LPL, this must deal with invalidated processes
   */
  rebase = (old_process: Process, new_process: Process) => {
    // Same as resolve but only up to old_process
    const fpl: FPL[] = [this]
    let head: FPL | undefined = this
    while (head.parent !== undefined) {
      if (head.process === old_process) break
      head = head.parent
      fpl.push(head)
    }
    // if we didn't end up on the old_proces, it's not in the FPL
    if (head.process !== old_process) {
      throw new Error(`${old_process} not found in FPL`)
    }
    // Replace old_process with new process
    const update = { [head.process.id]: new_process }
    // Walk forward updating processes, replacing any dependencies with the updated one
    head = new FPL(new_process, head.parent)
    const rebased = head
    fpl.pop()
    fpl.reverse()
    for (const el of fpl) {
      const new_proc = new Process(
        (update[el.process.id] || el.process).type,
        (update[el.process.id] || el.process).data,
        dict.init(
          dict.items(el.process.inputs).map(({ key, value }) =>
            ({ key, value: update[value.id] || value })
          )
        ),
        el.process.db,
      )
      head = new FPL(new_proc, head)
      update[el.process.id] = new_proc
    }
    return { rebased, head }
  }
}

type DatabaseKeyedTables = { process: Process, fpl: FPL, data: Data, resolved: Resolved }
type DatabaseListenCallback = <T extends keyof DatabaseKeyedTables>(table: T, record: DatabaseKeyedTables[T]) => void

export const IdOrDataC = z.union([
  z.object({ id: z.string() }),
  z.object({ type: z.string(), value: z.string() }),
])
export type IdOrData = z.infer<typeof IdOrDataC>

export type IdOrProcess = { id: string }
| { type: string, inputs: Record<string, IdOrProcess> }
| { type: string, data: IdOrData, inputs: Record<string, IdOrProcess> }

export const IdOrProcessC: z.ZodType<IdOrProcess> = z.lazy(() => z.union([
  z.object({
    id: z.string(),
  }),
  z.object({
    type: z.string(),
    data: IdOrDataC,
    inputs: z.record(z.string(), IdOrProcessC),
  }),
  z.object({
    type: z.string(),
    inputs: z.record(z.string(), IdOrProcessC),
  }),
]))

/**
 * The database maintains records for each table type
 */
export class Database {
  private processTable: Record<string, Process> = {}
  private fplTable: Record<string, FPL> = {}
  private resolvedTable: Record<string, Resolved> = {}
  private dataTable: Record<string, Data> = {}
  private listeners: Record<number, DatabaseListenCallback> = {}
  private id = 0

  /**
   * Listen to changes on the DB
   */
  listen = (cb: DatabaseListenCallback) => {
    const id = this.id++
    this.listeners[id] = cb
    return () => delete this.listeners[id]
  }
  /**
   * Inform listeners of changes to the DB
   */
  private notify = <T extends keyof DatabaseKeyedTables>(table: T, record: DatabaseKeyedTables[T]) => {
    for (const listener of Object.values(this.listeners)) {
      listener(table, record)
    }
  }

  resolveProcess = (process: IdOrProcess): Process => {
    if ('id' in process) {
      return this.getProcess(process.id) as Process
    } else {
      return this.upsertProcess(new Process(
        process.type,
        'data' in process ? process.data !== undefined ? this.resolveData(process.data) : undefined : undefined,
        dict.init(dict.items(process.inputs).map(({ key, value }) => ({ key, value: this.resolveProcess(value) })))
      ))
    }
  }
  getProcess = (id: string) => {
    return this.processTable[id] as Process | undefined
  }
  upsertProcess = (process: Process) => {
    if (!(process.id in this.processTable)) {
      process.db = this
      if (process.data !== undefined) {
        process.data = this.upsertData(process.data)
      }
      this.processTable[process.id] = process
      this.notify('process', process)
    }
    return this.processTable[process.id] as Process
  }

  resolveData = (data: IdOrData) => {
    if ('id' in data) {
      return this.getData(data.id) as Data
    } else {
      return this.upsertData(new Data(data.type, data.value))
    }
  }
  getData = (id: string) => {
    return this.dataTable[id] as Data | undefined
  }
  upsertData = (data: Data) => {
    if (!(data.id in this.dataTable)) {
      this.dataTable[data.id] = data
      this.notify('data', data)
    }
    return this.dataTable[data.id] as Data
  }

  getResolved = (id: string) => {
    return this.resolvedTable[id] as Resolved | undefined
  }
  /**
   * Like getResolved but if it's not in the db yet, wait for it
   */
  awaitResolved = async (id: string) => {
    if (!(id in this.resolvedTable)) {
      await new Promise<void>((resolve, reject) => {
        const unsub = this.listen((table, record) => {
          if (table === 'resolved' && record.id === id) {
            resolve()
            unsub()
          }
        })
      })
    }
    return this.resolvedTable[id] as Resolved
  }
  upsertResolved = (resolved: Resolved) => {
    if (!(resolved.id in this.resolvedTable)) {
      if (resolved.data) {
        this.upsertData(resolved.data)
      }
      this.resolvedTable[resolved.id] = resolved
      this.notify('resolved', resolved)
    }
    return this.resolvedTable[resolved.id] as Resolved
  }

  getFPL = (id: string) => {
    return this.fplTable[id] as FPL | undefined
  }
  upsertFPL = (fpl: FPL) => {
    if (!(fpl.id in this.fplTable)) {
      this.fplTable[fpl.id] = fpl
      fpl.process = this.upsertProcess(fpl.process)
      if (fpl.parent !== undefined) {
        fpl.parent = this.upsertFPL(fpl.parent)
      }
      this.notify('fpl', fpl)
    }
    return this.fplTable[fpl.id] as FPL
  }

  dump = () => {
    return {
      dataTable: this.dataTable,
      processTable: this.processTable,
      resolvedTable: this.resolvedTable,
      fplTable: this.fplTable,
    }
  }
}
