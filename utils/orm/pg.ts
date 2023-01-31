import { Codec, Decoded } from '@/spec/codec'
import { TypedSchema } from '@/spec/sql'
import * as pg from 'pg'
import * as dict from '@/utils/dict'
import { Create, DbTable, Delete, Find, FindMany, Update, Where } from './common'

/**
 * Easy prepared statement building.
 * Usage:
 * prepare(subst => `select * from item where id = ${subst('whatever')}`) == {
 *   statement: select * from item where id = $1
 *   vars: ['whatever']
 * }
 */
function prepare(statement_builder: (subst: (value: any) => string) => string) {
  let i = 1
  const vars: any[] = []
  const subst = (value: any) => {
    vars.push(value)
    return `$${i++}`
  }
  const statement = statement_builder(subst)
  return { statement, vars }
}

export class PgTable<T extends { id: Codec<string, string> }> implements DbTable<T> {
  constructor(public table: TypedSchema<T>, private pool: pg.Pool) {}

  private raw = async (statement_builder: (subst: (value: any) => string) => string) => {
    const { statement, vars } = prepare(statement_builder)
    return await this.pool.query(statement, vars)
  }
  
  create = async (create: Create<T>)=> {
    const columns = Object.keys(this.table.field_codecs).filter(col => col in create.data) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing values')
    const results = await this.raw(subst => `
      insert into ${JSON.stringify(this.table.name)} (${columns.map(col => JSON.stringify(col)).join(", ")})
      select * from
      jsonb_to_recordset(${subst(JSON.stringify([
        dict.init(columns.map((key) => ({ key, value: this.table.field_codecs[key].encode(create.data[key] as Decoded<T[keyof T]>) })))
      ]))}::jsonb) as t(
        ${columns.map(col => `${JSON.stringify(col)} ${this.table.field_sql[col]}`).join(',')}
      )
      ${this.table.extra_insert_sql || ''}
      returning *;
    `)
    if (results.rowCount !== 1) throw new Error('Expected one got none')
    return this.table.codec.decode(results.rows[0])
  }

  findUnique = async (find: Find<T>) => {
    const columns = Object.keys(this.table.field_codecs).filter(col => col in find.where) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing where')
    const results = await this.raw(subst => `
      select *
      from ${JSON.stringify(this.table.name)}
      where ${columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(find.where[key] as Decoded<T[keyof T]>))}`).join(' and ')};
    `)
    if (results.rowCount === 0) return null
    else if (results.rowCount > 1) throw new Error('Expected one got several')
    else return this.table.codec.decode(results.rows[0])
  }
  findMany = async (find: FindMany<T> = {}) => {
    const where: Where<T> = (find.where !== undefined) ? find.where : {}
    const columns = Object.keys(this.table.field_codecs).filter(col => col in where) as Array<keyof T>
    const results = await this.raw(subst => `
      select *
      from ${JSON.stringify(this.table.name)}
      ${columns.length > 0 ? columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(where[key] as Decoded<T[keyof T]>))}`).join(' and ') : ''};
    `)
    return results.rows.map(row => this.table.codec.decode(row))
  }
  update = async (update: Update<T>) => {
    const set_columns = Object.keys(this.table.field_codecs).filter(col => col in update.data) as Array<keyof T>
    if (set_columns.length <= 0) throw new Error('Missing data values')
    const where_columns = Object.keys(this.table.field_codecs).filter(col => col in update.where) as Array<keyof T>
    if (where_columns.length <= 0) throw new Error('Missing where values')
    const results = await this.raw(subst => `
      with rows as (
        update ${JSON.stringify(this.table.name)} as tbl
        set ${set_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode((update.data as any)[key]))}`).join(', ')}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(update.where[key] as Decoded<T[keyof T]>))}`).join(' and ')}
        returning tbl.*
      ) select one_and_only(*) from rows;
    `)
    if (results.rowCount === 0) return null
    else if (results.rowCount > 1) throw new Error('Unexpected output')
    else return this.table.codec.decode(results.rows[0])
  }
  updateMany = async (update: Update<T>) => {
    const set_columns = Object.keys(this.table.field_codecs).filter(col => col in update.data) as Array<keyof T>
    if (set_columns.length <= 0) throw new Error('Missing data values')
    const where_columns = Object.keys(this.table.field_codecs).filter(col => col in update.where) as Array<keyof T>
    if (where_columns.length <= 0) throw new Error('Missing where values')
    const results = await this.raw(subst => `
      with rows as (
        update ${JSON.stringify(this.table.name)}
        set ${set_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode((update.data as any)[key]))}`).join(', ')}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(update.where[key] as Decoded<T[keyof T]>))}`).join(' and ')};
        returning 1
      ) select count(*) as count from rows;
    `)
    return results.rows[0] as { count: number }
  }
  delete = async (delete_: Delete<T>) => {
    const where_columns = Object.keys(this.table.field_codecs).filter(col => col in delete_.where) as Array<keyof T>
    if (where_columns.length <= 0) throw new Error('Missing values')
    const results = await this.raw(subst => `
      with rows as (
        delete from ${JSON.stringify(this.table.name)}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(delete_.where[key] as Decoded<T[keyof T]>))}`).join(' and ')};
        returning *
      ) select one_and_only(*) from rows;
    `)
    if (results.rowCount === 0) return null
    else if (results.rowCount > 1) throw new Error('Unexpected output')
    else return this.table.codec.decode(results.rows[0])
  }
  deleteMany = async (delete_: Delete<T>) => {
    const where_columns = Object.keys(this.table.field_codecs).filter(col => col in delete_.where) as Array<keyof T>
    if (where_columns.length <= 0) throw new Error('Missing values')
    const results = await this.raw(subst => `
      with rows as (
        delete from ${JSON.stringify(this.table.name)}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(delete_.where[key] as Decoded<T[keyof T]>))}`).join(' and ')};
        returning 1
      ) select count(*) as count from rows;
    `)
    return results.rows[0] as { count: number }
  }
}
