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
import { z } from 'zod'
import { Db } from '@/utils/orm'
import * as fpprgSchema from '@/db/fpprg'
import { TypedSchema } from '@/spec/sql'

/**
 * This timeout error is used to ensure we don't wait too long for dependencies
 *  fortunately, even if it occurs the job will requeue still making progress.
 */
export class TimeoutError extends Error {}

/**
 * A process is unique by its own data (configuration) and process outputs
 *  which come before it. In that way, Processes operate like their own
 *  DAG alongside the FPL.
 */
export class Process {
  id: string
  resolved: Resolved | undefined = undefined

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
    /**
     * Whether or not this has been persisted to the db
     */
    public persisted = false,
  ) {
    this.id = uuid([type, data, dict.items(inputs).map(({ key, value }) => ({ key, value: value.id }))])
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
    return dict.init(
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

  constructor(public process: Process, public data: Data | undefined, public persisted = false) {
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
  constructor(public process: Process, public parent: FPL | undefined = undefined, public persisted = false) {
    this.id = uuid([process.id, parent ? parent.id : null])
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

/**
 * Content-addressable storage
 */
export class Data {
  id: string

  constructor(public type: string, public value: string, public persisted = false) {
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

export type DatabaseKeyedTables = { process: Process, fpl: FPL, data: Data, resolved: Resolved }
export type DatabaseListenCallback = <T extends keyof DatabaseKeyedTables>(table: T, record: DatabaseKeyedTables[T]) => void

type Schema = Pick<typeof fpprgSchema, 'data' | 'process' | 'process_input' | 'resolved' | 'process_complete' | 'fpl'>
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

  constructor(public db: Db<DbSchemaT>) {}

  resolveProcess = async (process: IdOrProcess): Promise<Process> => {
    if ('id' in process) {
      return await this.getProcess(process.id) as Process
    } else {
      return await this.upsertProcess(new Process(
        process.type,
        'data' in process ? process.data !== undefined ? await this.resolveData(process.data) : undefined : undefined,
        dict.init(await Promise.all(dict.items(process.inputs).map(async ({ key, value }) => ({ key, value: await this.resolveProcess(value) }))))
      ))
    }
  }
  getProcess = async (id: string, force = false) => {
    if (!(id in this.processTable) || force) {
      const result = await this.db.objects.process_complete.findUnique({ where: { id } })
      if (result !== null) {
        const process = new Process(
          result.type,
          result.data !== null ? await this.getData(result.data) : undefined,
          dict.init(await Promise.all(dict.items(result.inputs).map(async ({ key, value }) => ({ key, value: (await this.getProcess(value)) as Process })))),
          this,
          true,
        )
        this.processTable[id] = process
        // this process was already resolved?
        // save that information, saves a lookup or two
        if (result.resolved) {
          if (!(id in this.resolvedTable)) {
            const resolved = new Resolved(
              process,
              result.output !== null ? (await this.getData(result.output)) : undefined,
            )
            this.resolvedTable[id] = resolved
            process.resolved = resolved
          }
        }
      }
    }
    return this.processTable[id] as Process | undefined
  }
  upsertProcess = async (process: Process) => {
    if (!process.persisted) {
      if (!(process.id in this.processTable)) {
        process.db = this
        // TODO: get these things done in the context manager
        //       in a single transaction
        if (process.data !== undefined) {
          process.data = await this.upsertData(process.data)
        }
        for (const key in process.inputs) {
          const value = process.inputs[key]
          if (!value.persisted) {
            await this.upsertProcess(value)
          }
        }
        process.persisted = true
        this.processTable[process.id] = process
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
        //     inputs: dict.items(process.inputs).map(({ key, value }) => ({
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
      }
    }
    return this.processTable[process.id] as Process
  }

  resolveData = async (data: IdOrData) => {
    if ('id' in data) {
      return await this.getData(data.id) as Data
    } else {
      return await this.upsertData(new Data(data.type, data.value))
    }
  }
  getData = async (id: string) => {
    if (!(id in this.dataTable)) {
      const result = await this.db.objects.data.findUnique({ where: { id } })
      if (result !== null) {
        this.dataTable[result.id] = new Data(
          result.type,
          result.value,
          true,
        )
      }
    }
    return this.dataTable[id] as Data | undefined
  }
  upsertData = async (data: Data) => {
    if (!data.persisted) {
      if (!(data.id in this.dataTable)) {
        data.persisted = true
        this.dataTable[data.id] = data
        await this.db.objects.data.upsert({
          where: { id: data.id },
          create: {
            id: data.id,
            type: data.type,
            value: data.value,
          }
        })
      }
    }
    return this.dataTable[data.id] as Data
  }

  getResolved = async (id: string) => {
    if (!(id in this.resolvedTable)) {
      const result = await this.db.objects.resolved.findUnique({ where: { id } })
      if (result !== null) {
        const process = (await this.getProcess(result.id)) as Process
        const resolved = new Resolved(
          process,
          result.data !== null ? (await this.getData(result.data)) as Data : undefined,
          true,
        )
        process.resolved = resolved
        this.resolvedTable[result.id] = resolved
      }
    }
    return this.resolvedTable[id] as Resolved | undefined
  }
  /**
   * Like getResolved but if it's not in the db yet, request & wait for it
   */
  awaitResolved = async (id: string) => {
    const resolved = await this.getResolved(id)
    if (!resolved) {
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
        console.debug(`requesting ${id}`)
        await this.db.send('work-queue', { id })
        setTimeout(() => {
          if (!ctx.resolved) {
            ctx.resolved = true
            unsub()
            reject(new TimeoutError())
          }
        }, 5000)
      })
    }
    return this.resolvedTable[id] as Resolved
  }
  upsertResolved = async (resolved: Resolved) => {
    if (!resolved.persisted) {
      if (!(resolved.id in this.resolvedTable)) {
        if (resolved.data) {
          await this.upsertData(resolved.data)
        }
        resolved.persisted = true
        this.resolvedTable[resolved.id] = resolved
        resolved.process.resolved = resolved
        await this.db.objects.resolved.upsert({
          where: { id: resolved.id },
          create: {
            id: resolved.id,
            data: resolved.data !== undefined ? resolved.data.id : null,
          },
        })
      }
    }
    return this.resolvedTable[resolved.id] as Resolved
  }

  getFPL = async (id: string) => {
    if (!(id in this.fplTable)) {
      const result = await this.db.objects.fpl.findUnique({ where: { id } })
      if (result !== null) {
        this.fplTable[result.id] = new FPL(
          (await this.getProcess(result.process)) as Process,
          result.parent !== null ? (await this.getFPL(result.parent)) as FPL : undefined,
          true,
        )
      }
    }
    return this.fplTable[id] as FPL | undefined
  }
  upsertFPL = async (fpl: FPL) => {
    if (!fpl.persisted) {
      if (!(fpl.id in this.fplTable)) {
        fpl.process = await this.upsertProcess(fpl.process)
        if (fpl.parent !== undefined) {
          fpl.parent = await this.upsertFPL(fpl.parent)
        }
        fpl.persisted = true
        this.fplTable[fpl.id] = fpl
        await this.db.objects.fpl.upsert({
          where: { id: fpl.id },
          create: {
            id: fpl.id,
            process: fpl.process.id,
            parent: fpl.parent !== undefined ? fpl.parent.id : null,
          }
        })
      }
    }
    return this.fplTable[fpl.id] as FPL
  }

  dump = async () => {
    return {
      dataTable: this.dataTable,
      processTable: this.processTable,
      resolvedTable: this.resolvedTable,
      fplTable: this.fplTable,
    }
  }
}
