import type { TypedSchema, TypedSchemaRecord } from '@/spec/sql'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { Find, Create, Delete, Update, DbTable, FindMany, Upsert, DbDatabase } from './common'

export class MemoryDatabase implements DbDatabase {
  public objects: any
  private listeners: Record<number, (evt: string, data: unknown) => void> = {}
  private id = 0

  constructor(public data: Record<string, Record<string, Record<string, unknown>>> = {}) {}

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

export class MemoryTable<T extends {}> implements DbTable<T> {
  constructor(public table: TypedSchema<T>, private db: MemoryDatabase) {
    if (!(this.table.name in this.db.data)) this.db.data[this.table.name] = {}
  }

  private ensureCreate = async () => {
    if (this.table.js !== undefined) throw new Error('Cannot create on view')
  }
  private ensureFind = async () => {
    if (this.table.js !== undefined) {
      for (const data of (await this.table.js(this.db))) {
        const pk = this.table.field_pk.map((key) => `${data[key]}`).join('-')
        this.db.data[this.table.name][pk] = data
      }
    }
  }
  private ensureMutate = async () => {
    if (this.table.js !== undefined) throw new Error('Cannot modify view')
  }

  create = async (create: Create<T>)=> {
    await this.ensureCreate()
    const data = {...create.data}
    dict.items(this.table.field_default)
      .forEach(({ key, value: field_default }) => {
        if (typeof field_default === 'function' && (!(key in data) || data[key] === undefined || data[key] === null)) {
          data[key] = field_default()
        }
      })
    const pk = this.table.field_pk.map((key) => `${data[key]}`).join('-')
    this.db.data[this.table.name][pk] = data
    await this.db.notify(`insert:${this.table.name}`, { id: pk })
    return this.db.data[this.table.name][pk] as TypedSchemaRecord<TypedSchema<T>>
  }
  findUnique = async (find: Find<T>) => {
    await this.ensureFind()
    if (array.all(this.table.field_pk.map(key => key in find.where && find.where[key]))) {
      const pk = this.table.field_pk.map(key => `${find.where[key]}`).join('-')
      return (this.db.data[this.table.name][pk] || null) as TypedSchemaRecord<TypedSchema<T>> | null
    } else {
      const result = Object.values(this.db.data[this.table.name]).filter(record =>
        array.all(dict.items(find.where).map(({ key, value }) => record[key as string] === value))
      )
      if (result.length === 0) return null
      else if (result.length !== 1) throw new Error('Not unique')
      return result[0] as TypedSchemaRecord<TypedSchema<T>> | null
    }
  }
  findMany = async (find: FindMany<T> = {}) => {
    await this.ensureFind()
    return Object.values(this.db.data[this.table.name]).filter(record =>
      find.where ? array.all(dict.items(find.where).map(({ key, value }) => record[key as string] === value)) : true
    ) as Array<TypedSchemaRecord<TypedSchema<T>>>
  }
  update = async (update: Update<T>) => {
    await this.ensureMutate()
    const record = await this.findUnique({ where: update.where })
    if (record === null) return null
    Object.assign(record, update.data)
    return record as TypedSchemaRecord<TypedSchema<T>>
  }
  upsert = async (upsert: Upsert<T>) => {
    await this.ensureMutate()
    const record = upsert.update ? await this.update({ where: upsert.where, data: upsert.update }) : null
    if (record === null) {
      return await this.create({ data: upsert.create })
    } else {
      return record
    }
  }
  delete = async (delete_: Delete<T>) => {
    await this.ensureMutate()
    const record = await this.findUnique({ where: delete_.where })
    if (record === null) return null
    const pk = this.table.field_pk.map(key => `${record[key]}`).join('-')
    delete this.db.data[this.table.name][pk]
    return record as TypedSchemaRecord<TypedSchema<T>>
  }
}
