import type { Codec } from '@/spec/codec'
import type { TypedSchema, TypedSchemaRecord } from '@/spec/sql'
import { v4 as uuidv4 } from 'uuid'
import * as dict from '@/utils/dict'
import { Find, Create, Delete, Update, DbTable, FindMany } from './common'

const any = <T extends {}>(L: T[]) => L.some(el => el)
const all = <T extends {}>(L: T[]) => !L.some(el => !el)

export class MemoryTable<T extends { id: Codec<string, string> }> implements DbTable<T> {
  constructor(public table: TypedSchema<T>, private data: Record<string, Record<string | number | symbol, Record<string | number | symbol, unknown > & { id: string }>> = {}) {
    if (!(this.table.name in this.data)) this.data[this.table.name] = {}
  }

  create = async (create: Create<T>)=> {
    const id = create.data.id || uuidv4()
    this.data[this.table.name][id] = { ...create.data, id }
    return this.data[this.table.name][id] as TypedSchemaRecord<TypedSchema<T>>
  }
  findUnique = async (find: Find<T>) => {
    if (Object.keys(find.where).length === 1 && find.where.id) {
      return (this.data[this.table.name][find.where.id] || null) as TypedSchemaRecord<TypedSchema<T>> | null
    } else {
      const result = Object.values(this.data[this.table.name]).filter(record =>
        all(dict.items(find.where).map(({ key, value }) => record[key] === value))
      )
      if (result.length === 0) return null
      else if (result.length !== 1) throw new Error('Not unique')
      return result[0] as TypedSchemaRecord<TypedSchema<T>> | null
    }
  }
  findMany = async (find: FindMany<T> = {}) => {
    return Object.values(this.data[this.table.name]).filter(record =>
      find.where ? all(dict.items(find.where).map(({ key, value }) => record[key] === value)) : true
    ) as Array<TypedSchemaRecord<TypedSchema<T>>>
  }
  update = async (update: Update<T>) => {
    const record = await this.findUnique({ where: update.where })
    if (record === null) return null
    Object.assign(record, update.data)
    return record as TypedSchemaRecord<TypedSchema<T>>
  }
  delete = async (delete_: Delete<T>) => {
    const record = await this.findUnique({ where: delete_.where })
    if (record === null) return null
    delete this.data[this.table.name][record.id as string]
    return record as TypedSchemaRecord<TypedSchema<T>>
  }
}
