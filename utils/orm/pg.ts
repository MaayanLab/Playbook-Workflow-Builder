import { Decoded } from '@/spec/codec'
import { TypedSchema } from '@/spec/sql'
import * as pg from 'pg'
import * as dict from '@/utils/dict'
import { Create, DbDatabase, DbTable, Delete, Find, FindMany, OrderBy, Update, Upsert, WhereMany } from './common'
import PgBoss from 'pg-boss'
import createSubscriber, { Subscriber } from 'pg-listen'
import * as db from '@/db/orm'
import { z } from 'zod'

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

export class PgDatabase implements DbDatabase {
  public objects: any
  public pool: pg.Pool
  private listeners: Record<string, Record<number, (data: unknown) => void>> = {}
  private id = 0
  private subscriber: Subscriber
  private boss: PgBoss

  constructor(connectionString: string) {
    this.pool = new pg.Pool({ connectionString: connectionString })
    this.subscriber = createSubscriber({ connectionString })
    this.boss = new PgBoss(connectionString)
    this.boss.on('error', error => console.error(error))
    this.subscriber.events.on('error', (error) => {
      console.error("Fatal database connection error:", error)
      process.exit(1)
    })
    this.boss.on('error', (error) => {
      console.error("PgBoss error:", error)
    })
    process.on("exit", () => {
      this.subscriber.close()
    })
    ;(async (self) => {
      await self.boss.start()
      await self.subscriber.connect()
      this.listen('distributed:on_insert', async (rawPayload) => {
        const payload = db.notify_insertion_trigger_payload.codec.decode(rawPayload)
        const { table, operation, id } = payload
        if (operation === 'INSERT') {
          this.notify(`insert:${table}`, { id })
        }
      })
    })(this).catch((error) => {
      console.error('Failed to initialize subscriber', error)
      process.exit(1)
    })
  }

  raw = async (statement_builder: (subst: (value: any) => string) => string) => {
    const { statement, vars } = prepare(statement_builder)
    return await this.pool.query(statement, vars)
  }

  listen = <T>(evt: string, cb: (data: T) => void) => {
    const id = this.id++
    if (!(evt in this.listeners)) {
      this.listeners[evt] = {}
      if (evt.startsWith('distributed:')) {
        const dbEvt = evt.slice(evt.indexOf(':')+1)
        this.subscriber.listenTo(dbEvt)
        this.subscriber.notifications.on(dbEvt, (rawPayload) => {
          this._notify(evt, rawPayload)
        })
      }
    }
    this.listeners[evt][id] = cb as (data: unknown) => void
    return () => {
      delete this.listeners[evt][id]
      if (dict.isEmpty(this.listeners[evt])) {
        delete this.listeners[evt]
        if (evt.startsWith('distributed:')) {
          const dbEvt = evt.slice(evt.indexOf(':'))
          this.subscriber.unlisten(dbEvt)
        }
      }
    }
  }
  /*
   * Send events with distributed: to the db
   * but don't send receive them here
   */
  notify = <T>(evt: string, data: T) => {
    if (evt.startsWith('distributed:')) {
      const dbEvt = evt.slice(evt.indexOf(':')+1)
      this.subscriber.notify(dbEvt, data)
    } else {
      this._notify(evt, data)
    }
  }

  _notify = <T>(evt: string, data: T) => {
    for (const listener of dict.values(this.listeners[evt] ?? {})) {
      listener(data)
    }
  }

  send = async (queue: string, work: { id: string, priority?: number }) => {
    await this.boss.send(queue, { id: work.id }, { singletonKey: work.id, priority: work.priority ?? 0 })
  }

  work = async (queue: string, opts: unknown, cb: (work: { data: { id: string } }) => Promise<void>) => {
    await this.boss.work(queue, opts as any, cb)
    return () => {
      this.boss.stop().catch(error => console.error(error))
    }
  }

  jobs = async () => {
    const results = await this.raw(subst => `
      select
        name,
        priority,
        data,
        state,
        createdon,
        completedon
      from pgboss.job
      where data is not null
    `)
    return z.array(z.object({
      name: z.string(),
      priority: z.number(),
      data: z.any(),
      state: z.string(),
      createdon: z.date(),
      completedon: z.date().nullable(),
    })).parse(results.rows)
  }
}

export class PgTable<T extends {}> implements DbTable<T> {
  constructor(public table: TypedSchema<T>, private db: PgDatabase) {}

  create = async (create: Create<T>)=> {
    const columns = dict.keys(this.table.field_codecs).filter(col => col in create.data) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing values')
    const results = await this.db.raw(subst => `
      insert into ${JSON.stringify(this.table.name)} (${columns.map(col => JSON.stringify(col)).join(", ")})
      select * from
      jsonb_to_recordset(${subst(JSON.stringify([
        dict.init(columns.map((key) => ({ key, value: this.table.field_codecs[key].encode(create.data[key] as Decoded<T[keyof T]>) })))
      ]))}::jsonb) as t(
        ${columns.map(col => `${JSON.stringify(col)} ${this.table.field_sql[col]}`).join(',')}
      )
      returning *;
    `)
    if (results.rowCount !== 1) throw new Error('Expected one got none')
    return this.table.codec.decode(results.rows[0])
  }

  findUnique = async (find: Find<T>) => {
    const columns = dict.keys(this.table.field_codecs).filter(col => col in find.where) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing where')
    const results = await this.db.raw(subst => `
      select *
      from ${JSON.stringify(this.table.name)}
      where ${columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(find.where[key] as Decoded<T[keyof T]>))}`).join(' and ')};
    `)
    if (results.rowCount === 0 || results.rowCount === null) return null
    else if (results.rowCount > 1) throw new Error('Expected one got several')
    else return this.table.codec.decode(results.rows[0])
  }
  findMany = async (find: FindMany<T> = {}) => {
    const where: WhereMany<T> = (find.where !== undefined) ? find.where : {}
    const where_columns = dict.keys(this.table.field_codecs).filter(col => col in where) as Array<keyof T>
    const orderBy: OrderBy<T> = (find.orderBy !== undefined) ? find.orderBy : {}
    const orderBy_columns = dict.keys(this.table.field_codecs).filter(col => col in orderBy) as Array<keyof T>
    if (find.skip !== undefined && typeof find.skip !== 'number') throw new Error('Expected number for skip')
    if (find.take !== undefined && typeof find.take !== 'number') throw new Error('Expected number for take')
    const results = await this.db.raw(subst => `
      select *
      from ${JSON.stringify(this.table.name)}
      ${where_columns.length > 0 ?
        `where ${where_columns
          .map(key => {
            const value = where[key]
            if (typeof value === 'object' && value !== null && 'in' in value) {
              return `${JSON.stringify(key)} = any(${subst((value.in as Decoded<T[keyof T]>[]).map(this.table.field_codecs[key].encode))})`
            } else {
              return `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(value as Decoded<T[keyof T]>))}`
            }
          })
          .join(' and ')
        }`
        : ''}
      ${orderBy_columns.length > 0 ?
        `order by ${orderBy_columns
          .map(key => `${JSON.stringify(key)} ${orderBy[key] === 'desc' ? 'desc' : 'asc'}`)
          .join(', ')
        }`
        : ''}
      ${find.take ? `limit ${find.take}` : ''}
      ${find.skip ? `offset ${find.skip}` : ''}
      ;
    `)
    return results.rows.map(row => this.table.codec.decode(row))
  }
  update = async (update: Update<T>) => {
    const set_columns = dict.keys(this.table.field_codecs).filter(col => col in update.data) as Array<keyof T>
    if (set_columns.length <= 0) throw new Error('Missing data values')
    const where_columns = dict.keys(this.table.field_codecs).filter(col => col in update.where) as Array<keyof T>
    if (where_columns.length <= 0) throw new Error('Missing where values')
    const results = await this.db.raw(subst => `
      with rows as (
        update ${JSON.stringify(this.table.name)} as tbl
        set ${[
          ...set_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode((update.data as any)[key]))}`),
          ...this.table.on_set_extra.map(([key, value]) => `${key} = ${value}`),
        ].join(', ')}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(update.where[key] as Decoded<T[keyof T]>))}`).join(' and ')}
        returning tbl.*
      ), assertion as (
        select assert_is_one(count(*))
        from rows
      ) select * from rows;
    `)
    if (results.rowCount === 0 || results.rowCount === null) return null
    else if (results.rowCount > 1) throw new Error('Unexpected output')
    else return this.table.codec.decode(results.rows[0])
  }
  updateMany = async (update: Update<T>) => {
    const set_columns = dict.keys(this.table.field_codecs).filter(col => col in update.data) as Array<keyof T>
    if (set_columns.length <= 0) throw new Error('Missing data values')
    const where_columns = dict.keys(this.table.field_codecs).filter(col => col in update.where) as Array<keyof T>
    if (where_columns.length <= 0) throw new Error('Missing where values')
    const results = await this.db.raw(subst => `
      with rows as (
        update ${JSON.stringify(this.table.name)}
        set ${[
          ...set_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode((update.data as any)[key]))}`),
        ...this.table.on_set_extra.map(([key, value]) => `${key} = ${value}`),
        ].join(', ')}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(update.where[key] as Decoded<T[keyof T]>))}`).join(' and ')}
        returning 1
      ) select count(*) as count from rows;
    `)
    return results.rows[0] as { count: number }
  }
  upsert = async (upsert: Upsert<T>) => {
    const insert_columns = dict.keys(this.table.field_codecs).filter(col => col in upsert.create) as Array<keyof T>
    if (insert_columns.length <= 0) throw new Error('Missing values')
    const update_columns = (upsert.update ? dict.keys(this.table.field_codecs).filter(col => upsert.update && col in upsert.update) : []) as Array<keyof T>
    this.table.field_pk.forEach((key) => {
      if (!update_columns.includes(key as keyof T)) update_columns.push(key as keyof T)
    })
    const results = await this.db.raw(subst => `
      insert into ${JSON.stringify(this.table.name)} (${insert_columns.map(col => JSON.stringify(col)).join(", ")})
      select * from
      jsonb_to_recordset(${subst(JSON.stringify([
        dict.init(insert_columns.map((key) => ({ key, value: this.table.field_codecs[key].encode(upsert.create[key] as Decoded<T[keyof T]>) })))
      ]))}::jsonb) as t(
        ${insert_columns.map(col => `${JSON.stringify(col)} ${this.table.field_sql[col]}`).join(',')}
      )
      on conflict (${this.table.field_pk.map((key) => JSON.stringify(key)).join(',')})
      do update set ${[
        ...update_columns.map(col => `${JSON.stringify(col)} = EXCLUDED.${JSON.stringify(col)}`),
        ...this.table.on_set_extra.map(([key, value]) => `${key} = EXCLUDED.${value}`),
      ].join(',')}
      returning *;
    `)
    if (results.rowCount !== 1) throw new Error('Expected one got none')
    return this.table.codec.decode(results.rows[0])
  }
  delete = async (delete_: Delete<T>) => {
    const where_columns = dict.keys(this.table.field_codecs).filter(col => col in delete_.where) as Array<keyof T>
    if (where_columns.length <= 0) throw new Error('Missing values')
    const results = await this.db.raw(subst => `
      with rows as (
        delete from ${JSON.stringify(this.table.name)}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(delete_.where[key] as Decoded<T[keyof T]>))}`).join(' and ')}
        returning *
      ), assertion as (
        select assert_is_one(count(*))
        from rows
      ) select * from rows;
    `)
    if (results.rowCount === 0 || results.rowCount === null) return null
    else if (results.rowCount > 1) throw new Error('Unexpected output')
    else return this.table.codec.decode(results.rows[0])
  }
  deleteMany = async (delete_: Delete<T>) => {
    const where_columns = dict.keys(this.table.field_codecs).filter(col => col in delete_.where) as Array<keyof T>
    if (where_columns.length <= 0) throw new Error('Missing values')
    const results = await this.db.raw(subst => `
      with rows as (
        delete from ${JSON.stringify(this.table.name)}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(delete_.where[key] as Decoded<T[keyof T]>))}`).join(' and ')}
        returning 1
      ) select count(*) as count from rows;
    `)
    return results.rows[0] as { count: number }
  }
}
