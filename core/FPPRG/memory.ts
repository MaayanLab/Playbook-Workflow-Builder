/**
 * This is an in-memory database implementation of FPKRG
 */

import * as dict from '@/utils/dict'
import { Database, Process, FPL, Resolved, Data, DatabaseListenCallback, DatabaseKeyedTables, IdOrProcess, IdOrData } from '@/core/FPPRG/spec'

/**
 * The database maintains records for each table type
 */
export class MemoryDatabase implements Database {
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
    return () => {delete this.listeners[id]}
  }
  /**
   * Inform listeners of changes to the DB
   */
  private notify = async <T extends keyof DatabaseKeyedTables>(table: T, record: DatabaseKeyedTables[T]) => {
    for (const listener of Object.values(this.listeners)) {
      listener(table, record)
    }
  }

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
  getProcess = async (id: string) => {
    return this.processTable[id] as Process | undefined
  }
  upsertProcess = async (process: Process) => {
    if (!process.persisted) {
      if (!(process.id in this.processTable)) {
        process.db = this
        if (process.data !== undefined) {
          process.data = await this.upsertData(process.data)
        }
        process.persisted = true
        this.processTable[process.id] = process
        await this.notify('process', process)
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
    return this.dataTable[id] as Data | undefined
  }
  upsertData = async (data: Data) => {
    if (!data.persisted) {
      if (!(data.id in this.dataTable)) {
        data.persisted = true
        this.dataTable[data.id] = data
        await this.notify('data', data)
      }
    }
    return this.dataTable[data.id] as Data
  }

  getResolved = async (id: string) => {
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
  upsertResolved = async (resolved: Resolved) => {
    if (!resolved.persisted) {
      if (!(resolved.id in this.resolvedTable)) {
        if (resolved.data) {
          await this.upsertData(resolved.data)
        }
        resolved.persisted = true
        this.resolvedTable[resolved.id] = resolved
        await this.notify('resolved', resolved)
      }
    }
    return this.resolvedTable[resolved.id] as Resolved
  }

  getFPL = async (id: string) => {
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
        await this.notify('fpl', fpl)
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
