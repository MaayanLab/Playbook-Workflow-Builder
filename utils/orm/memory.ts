import type { Codec } from '@/spec/codec'
import type { TypedSchema, TypedSchemaRecord } from '@/spec/sql'
import { v4 as uuidv4 } from 'uuid'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { Find, Create, Delete, Update, DbTable, FindMany, Upsert, DbDatabase } from './common'

export class MemoryDatabase implements DbDatabase {
  private listeners: Record<number, (evt: string, data: unknown) => void> = {}
  private id = 0

  constructor(public data: Record<string, Record<string | number | symbol, Record<string | number | symbol, unknown > & { id: string }>> = {}) {}

  listen = (cb: (evt: string, data: unknown) => void) => {
    const id = this.id++
    this.listeners[id] = cb
    return () => {delete this.listeners[id]}
  }

  notify = async (evt: string, data: unknown) => {
    for (const listener of Object.values(this.listeners)) {
      listener(evt, data)
    }
  }

  send = async (queue: string, work: unknown) => {
    this.notify(`boss:${queue}`, work)
  }

  work = async (queue: string, opts: unknown, cb: (work: unknown) => Promise<void>) => {
    return this.listen(async (evt, work: unknown) => {
      if (evt === `boss:${queue}`) await cb({ data: work })
    })
  }
}

export class MemoryTable<T extends { id: Codec<string, string> }> implements DbTable<T> {
  constructor(public table: TypedSchema<T>, private db: MemoryDatabase) {
    if (!(this.table.name in this.db.data)) this.db.data[this.table.name] = {}
  }

  create = async (create: Create<T>)=> {
    const id = create.data.id || uuidv4()
    this.db.data[this.table.name][id] = { ...create.data, id }
    await this.db.notify(`insert:${this.table.name}`, { id })
    return this.db.data[this.table.name][id] as TypedSchemaRecord<TypedSchema<T>>
  }
  findUnique = async (find: Find<T>) => {
    if (Object.keys(find.where).length === 1 && find.where.id) {
      return (this.db.data[this.table.name][find.where.id] || null) as TypedSchemaRecord<TypedSchema<T>> | null
    } else {
      const result = Object.values(this.db.data[this.table.name]).filter(record =>
        array.all(dict.items(find.where).map(({ key, value }) => record[key] === value))
      )
      if (result.length === 0) return null
      else if (result.length !== 1) throw new Error('Not unique')
      return result[0] as TypedSchemaRecord<TypedSchema<T>> | null
    }
  }
  findMany = async (find: FindMany<T> = {}) => {
    return Object.values(this.db.data[this.table.name]).filter(record =>
      find.where ? array.all(dict.items(find.where).map(({ key, value }) => record[key] === value)) : true
    ) as Array<TypedSchemaRecord<TypedSchema<T>>>
  }
  update = async (update: Update<T>) => {
    const record = await this.findUnique({ where: update.where })
    if (record === null) return null
    Object.assign(record, update.data)
    return record as TypedSchemaRecord<TypedSchema<T>>
  }
  upsert = async (upsert: Upsert<T>) => {
    const record = upsert.update ? await this.update({ where: upsert.where, data: upsert.update }) : null
    if (record === null) {
      return await this.create({ data: upsert.create })
    } else {
      return record
    }
  }
  delete = async (delete_: Delete<T>) => {
    const record = await this.findUnique({ where: delete_.where })
    if (record === null) return null
    delete this.db.data[this.table.name][record.id as string]
    return record as TypedSchemaRecord<TypedSchema<T>>
  }
}
