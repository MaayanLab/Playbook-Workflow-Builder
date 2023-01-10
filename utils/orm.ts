import { TypedSchema } from '@/spec/sql'
import * as pg from 'pg'

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
  private schemas: Record<string, TypedSchema<unknown>> = {}
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
  insert = async <C extends Record<string, unknown>>(table: TypedSchema<C>, ...data: C[]) => {
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
        const columns = Object.keys(schema.field_types)
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
