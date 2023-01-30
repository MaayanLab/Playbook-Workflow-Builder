import { Codec, Decoded, Encoded } from '@/spec/codec'
import { TableSchema, TypedSchema } from '@/spec/sql'
import * as pg from 'pg'
import { z } from 'zod'
import * as dict from '@/utils/dict'

/**
 * Easy prepared statement building.
 * Usage:
 * prepare(subst => `select * from item where id = ${subst('whatever')}`) == {
 *   statement: select * from item where id = $1
 *   vars: ['whatever']
 * }
 */
export function prepare(statement_builder: (subst: (value: any) => string) => string) {
  let i = 1
  const vars: any[] = []
  const subst = (value: any) => {
    vars.push(value)
    return `$${i++}`
  }
  const statement = statement_builder(subst)
  return { statement, vars }
}

export class ResultHelper<D, E> {
  constructor(private raw: pg.QueryResult, private codec: Codec<D, E>) {}
  one = () => {
    if (this.raw.rowCount === 0)  throw new Error('Expected one, received none')
    else if (this.raw.rowCount === 1) return this.codec.decode(this.raw.rows[0])
    else throw new Error('Expected one, received more than one')
  }
  oneOrNone = () => {
    if (this.raw.rowCount === 0) return null
    else if (this.raw.rowCount === 1) return this.codec.decode(this.raw.rows[0])
    else throw new Error('Expected one, received more than one')
  }
  iterator = () => {
    return this.raw.rows.map(row => this.codec.decode(row))
  }
}

const identity_codec: Codec<any, any> = { encode: T => T, decode: T => T }
export class ClientHelper {
  constructor(private client: pg.Client | pg.Pool) {}
  query = async (statement_builder: (subst: (value: any) => string) => string, codec = identity_codec) => {
    const { statement, vars } = prepare(statement_builder)
    return new ResultHelper(await this.client.query(statement, vars), codec)
  }
  insert = async <T>(schema: TypedSchema<T>, record: Partial<{ [K in keyof T]: Decoded<T[K]> }>) => {
    const columns = Object.keys(schema.field_codecs).filter(col => col in record) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing values')
    return await this.query(subst => `
      insert into ${JSON.stringify(schema.name)} (${columns.map(col => JSON.stringify(col)).join(", ")})
      select * from
      jsonb_to_recordset(${subst(
        dict.init(columns.map((key) => ({ key, value: schema.field_codecs[key].encode(record[key] as Decoded<T[keyof T]>) })))
      )}) as t(
        ${columns.map(col => `${JSON.stringify(col)} ${schema.field_sql[col]}`).join(',')}
      )
      ${schema.extra_insert_sql || ''};
    `)
  }
  insertReturning = async <T>(schema: TypedSchema<T>, record: Partial<{ [K in keyof T]: Decoded<T[K]> }>) => {
    const columns = Object.keys(schema.field_codecs).filter(col => col in record) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing values')
    return await this.query(subst => `
      insert into ${JSON.stringify(schema.name)} (${columns.map(col => JSON.stringify(col)).join(", ")})
      select * from
      jsonb_to_recordset(${subst(
        dict.init(columns.map((key) => ({ key, value: schema.field_codecs[key].encode(record[key] as Decoded<T[keyof T]>) })))
      )}) as t(
        ${columns.map(col => `${JSON.stringify(col)} ${schema.field_sql[col]}`).join(',')}
      )
      ${schema.extra_insert_sql || ''}
      returning *;
    `, schema.codec)
  }
  updateReturning = async <T>(schema: TypedSchema<T>, values: { id: string } & Partial<Omit<{ [K in keyof T]: Decoded<T[K]> }, 'id'>>) => {
    const columns = Object.keys(schema.field_codecs).filter(col => col in values && col !== 'id') as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing values')
    return await this.query(subst => `
      update ${JSON.stringify(schema.name)}
      set ${columns.map(key => `${JSON.stringify(key)} = ${subst(schema.field_codecs[key].encode((values as any)[key]))}`).join(', ')}
      where "id" = ${subst(values.id)}
      returning *;
    `, schema.codec)
  }
  selectWhere = async <T>(schema: TypedSchema<T>, where: Partial<{ [K in keyof T]: Decoded<T[K]> }>) => {
    const columns = Object.keys(schema.field_codecs).filter(col => col in where) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing where')
    return await this.query(subst => `
      select *
      from ${JSON.stringify(schema.name)}
      where ${columns.map(key => `${JSON.stringify(key)} = ${subst(schema.field_codecs[key].encode(where[key] as Decoded<T[keyof T]>))}`).join(' and ')};
    `, schema.codec)
  }
  deleteWhere = async <T>(schema: TypedSchema<T>, where: Partial<{ [K in keyof T]: Decoded<T[K]> }>) => {
    const columns = Object.keys(schema.field_codecs).filter(col => col in where) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing where')
    return await this.query(subst => `
      delete
      from ${JSON.stringify(schema.name)}
      where ${columns.map(key => `${JSON.stringify(key)} = ${subst(schema.field_codecs[key].encode(where[key] as Decoded<T[keyof T]>))}`).join(' and ')};
    `)
  }
  deleteWhereReturning = async <T>(schema: TypedSchema<T>, where: Partial<{ [K in keyof T]: Decoded<T[K]> }>) => {
    const columns = Object.keys(schema.field_codecs).filter(col => col in where) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing where')
    return await this.query(subst => `
      delete
      from ${JSON.stringify(schema.name)}
      where ${columns.map(key => `${JSON.stringify(key)} = ${subst(schema.field_codecs[key].encode(where[key] as Decoded<T[keyof T]>))}`).join(' and ')}
      returning *;
    `, schema.codec)
  }
}

/**
 * BatchInsertion lets us insert type-safe records in the way we find easiest
 *  and will send many insertions at a time when enough records have accumulated
 *  or when we're done.
 * 
 * Usage:
 * const pool = pg.Pool()
 * 
 * // create schema for the table containing sql & zod schema
 * const x_tbl = Table.create('x_tbl').field('x', 'int', '', z.number()).build()
 * const y_tbl = Table.create('y_tbl').field('y', 'int', '', z.number()).build()
 * const z_tbl = Table.create('z_tbl').field('z', 'int', '', z.number()).build()
 * 
 * // perform insertions in BatchInsertion context
 * await BatchInsertion.context(pool, async (ctx) => {
 *   // we can insert records one at a time with different tables
 *   for (let i = 0; i < 10; i++) {
 *     await ctx.insert(x_tbl, { x: i });
 *     await ctx.insert(y_tbl, { y: i*i });
 *   }
 *   // we can insert multiple records in one line
 *   await ctx.insert(z_tbl, ...[1, 1, 3, 5].map(z => ({ z })))
 * })
 */
export class BatchInsertion {
  private schemas: Record<string, TableSchema<unknown>> = {}
  private pending: Record<string, Record<string, unknown>[]> = {}

  constructor(private client: pg.PoolClient, private queue_size: number = 100) {}

  /**
   * BatchInsertion full lifecycle
   */
  static async context(pool: pg.Pool, cb: (ctx: BatchInsertion) => Promise<void>, opts?: { queue_size: number }) {
    const client = await pool.connect()
    try {
      await client.query('BEGIN')
      const ctx = new BatchInsertion(client, opts?.queue_size || 100)
      await cb(ctx)
      await ctx.flush()
      await client.query('COMMIT')
    } catch (e) {
      await client.query('ROLLBACK')
      throw e
    } finally {
      client.release()
    }
  }

  /**
   * Save record(s) for a given schema
   */
  insert = async <C extends Record<string, unknown>>(table: TableSchema<C>, ...data: C[]) => {
    if (!(table.name in this.schemas)) {
      this.schemas[table.name] = table
      this.pending[table.name] = []
    }
    for (const d of data) {
      this.pending[table.name].push(d)
      if (this.pending[table.name].length > this.queue_size) {
        await this.flush()
      }
    }
  }

  /**
   * Send pending records
   */
  flush = async () => {
    for (const table in this.schemas) {
      const schema = this.schemas[table]
      const pending = this.pending[table]
      if (pending.length === 0) {
        continue
      } else {
        const columns = Object.keys(schema.field_codecs)
        const { statement, vars } = prepare(subst => `
          insert into ${schema.name} (
            ${columns.map(field => JSON.stringify(field)).join(',')}
          )
          select * from
          jsonb_to_recordset(${subst(pending)}) as t(
            ${columns.map(field => `${JSON.stringify(field)} ${schema.field_sql}`).join(',')}
          )
          ${schema.extra_insert_sql || ''};
        `)
        await this.client.query(statement, vars)
        this.pending[table] = []
      }
    }
  }
}
