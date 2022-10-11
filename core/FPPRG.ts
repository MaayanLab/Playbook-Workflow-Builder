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

import { Hash } from "fast-sha256"
import * as dict from '@/utils/dict'

/**
 * Obtain a unique hash for a arbitrary data
 */
function sha256(data: any) {
  const h = new Hash()
  h.update(Buffer.from(JSON.stringify(data)))
  return Buffer.from(h.digest()).toString('hex')
}

/**
 * Content-addressable storage
 */
export class Data<T = unknown> {
  id: string

  constructor(public type: string, public value: T) {
    this.id = sha256([type, value])
  }

  toJSON = () => {
    return {
      '@id': this.id,
      '@type': this.type,
      '@value': this.value,
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
    public type: string,
    public inputs: Record<string|number, Process> = {},
    public db: Database = undefined,
  ) {
    this.id = sha256([type, dict.items(inputs).map(({ key, value }) => ({ key, value: value.id }))])
  }

  toJSON = (): { '@id': string, '@type': string, 'inputs': object } => {
    return {
      '@id': this.id,
      '@type': this.type,
      'inputs': dict.init(
        dict.items(this.inputs)
          .map(({ key, value }) => ({ key, value: value.toJSON() }))
      ) as object,
    }
  }

  /**
   * The outputs of all input processes
   */
  async inputs__outputs<O = unknown>(): Promise<Record<string | number | symbol, Data<O>>> {
    return dict.init<Data<O>>(
      await Promise.all(
        dict.items(this.inputs).map(async ({ key, value }) =>
          ({ key, value: await value.output<O>() })
        )
      )
    )
  }
  /**
   * The output of this process
   */
  async output<T = unknown>() {
    if (this.db === undefined) throw new Error('Process not attached to a db')
    const resolved = await this.db.awaitResolved<T>(this.id)
    return resolved.data
  }
}

/**
 * Data resolved by a process
 */
export class Resolved<T = unknown> {
  public id: string

  constructor(public process: Process, public data: Data<T>) {
    this.id = process.id
  }

  toJSON = () => {
    return {
      '@id': this.id,
      'data': this.data.toJSON(),
    }
  }
}

/**
 * Fully Persistent List of processes to maintain ordering
 */
export class FPL {
  id: string
  constructor(public process: Process, public parent: FPL | undefined = undefined) {
    this.id = sha256([process.id, parent ? parent.id : undefined])
  }

  toJSON = () => {
    return {
      '@id': this.id,
      'process': this.process.toJSON(),
    }
  }

  /**
   * Resolve the complete list from the FPL data structure
   */
  resolve() {
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
   * Extend an LPL with a new process
   */
  extend(process: Process) {
    return new FPL(process, this)
  }
  /**
   * Rebase an LPL, this must deal with invalidated processes
   */
  rebase(old_process: Process, new_process: Process)  {
    // Same as resolve but only up to old_process
    const fpl: FPL[] = [this]
    let head: FPL = this
    while (head.parent !== undefined) {
      if (head.process === old_process) break
      head = head.parent
      fpl.push(head)
    }
    fpl.reverse()
    // if we didn't end up on the old_proces, it's not in the FPL
    if (head.process !== old_process) {
      throw new Error(`${old_process} not found in FPL`)
    }
    // Replace old_process with new process
    const update = { [head.process.id]: new_process }
    head = new FPL(new_process, head.parent)
    // Walk forward updating processes, replacing any dependencies with the updated one
    for (const el of fpl) {
      const new_proc = new Process(
        el.process.type,
        dict.init(
          dict.items(el.process.inputs).map(({ key, value }) =>
            ({ key, value: update[value.id] || value })
          )
        )
      )
      head = new FPL(new_proc, head)
      update[el.process.id] = new_proc
    }
    return head
  }
}

type DatabaseKeyedTables = { process: Process, fpl: FPL, data: Data, resolved: Resolved }
type DatabaseListenCallback = <T extends keyof DatabaseKeyedTables>(table: T, record: DatabaseKeyedTables[T]) => void

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
  listen(cb: DatabaseListenCallback) {
    const id = this.id++
    this.listeners[id] = cb
    return () => delete this.listeners[id]
  }
  /**
   * Inform listeners of changes to the DB
   */
  private notify<T extends keyof DatabaseKeyedTables>(table: T, record: DatabaseKeyedTables[T]) {
    for (const listener of Object.values(this.listeners)) {
      listener(table, record)
    }
  }

  getProcess(id: string) {
    return this.processTable[id] as Process | undefined
  }
  upsertProcess(process: Process) {
    if (!(process.id in this.processTable)) {
      process.db = this
      this.processTable[process.id] = process
      this.notify('process', process)
    }
    return this.processTable[process.id]
  }

  getData<T = unknown>(id: string) {
    return this.dataTable[id] as Data<T> | undefined
  }
  upsertData<T>(data: Data<T>) {
    if (!(data.id in this.dataTable)) {
      this.dataTable[data.id] = data
      this.notify('data', data)
    }
    return this.dataTable[data.id]
  }

  getResolved<T = unknown>(id: string) {
    return this.resolvedTable[id] as Resolved<T> | undefined
  }
  /**
   * Like getResolved but if it's not in the db yet, wait for it
   */
  async awaitResolved<T = unknown>(id: string) {
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
    return this.resolvedTable[id] as Resolved<T>
  }
  upsertResolved<T>(resolved: Resolved<T>) {
    if (!(resolved.id in this.resolvedTable)) {
      this.upsertData(resolved.data)
      this.resolvedTable[resolved.id] = resolved
      this.notify('resolved', resolved)
    }
    return this.resolvedTable[resolved.id]
  }

  getFPL(id: string) {
    return this.fplTable[id]
  }
  upsertFPL(fpl: FPL) {
    if (!(fpl.id in this.fplTable)) {
      this.fplTable[fpl.id] = fpl
      fpl.process = this.upsertProcess(fpl.process)
      if (fpl.parent !== undefined) {
        fpl.parent = this.upsertFPL(fpl.parent)
      }
      this.notify('fpl', fpl)
    }
    return this.fplTable[fpl.id]
  }

  dump() {
    return {
      dataTable: this.dataTable,
      processTable: this.processTable,
      resolvedTable: this.resolvedTable,
      fplTable: this.fplTable,
    }
  }
}
