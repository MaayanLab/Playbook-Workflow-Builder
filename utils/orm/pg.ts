import { Codec, Decoded } from '@/spec/codec'
import { TypedSchema } from '@/spec/sql'
import * as pg from 'pg'
import * as dict from '@/utils/dict'
import { Create, DbDatabase, DbTable, Delete, Find, FindMany, Update, Upsert, Where } from './common'
import PgBoss from 'pg-boss'
import createSubscriber, { Subscriber } from 'pg-listen'
import * as db from '@/db/fpprg'

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
  private listeners: Record<number, (evt: string, data: unknown) => void> = {}
  private id = 0
  private subscriber: Subscriber
  private pool: pg.Pool
  private boss: PgBoss

  constructor(connectionString: string) {
    this.pool = new pg.Pool({ connectionString: connectionString })
    this.subscriber = createSubscriber({ connectionString })
    this.boss = new PgBoss(connectionString)
    this.boss.on('error', error => console.error(error))
    this.subscriber.notifications.on('on_insert', async (rawPayload) => {
      const payload = db.notify_insertion_trigger_payload.codec.decode(rawPayload)
      const { table, operation, id } = payload
      console.debug(`received on_insert from ${table}: ${id}`)
      if (operation === 'INSERT') {
        this.notify(`insert:${table}`, { id })
      }
    })
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
      await self.subscriber.listenTo('on_insert')
      console.log('ready')
    })(this).catch((error) => {
      console.error('Failed to initialize subscriber', error)
      process.exit(1)
    })
  }

  raw = async (statement_builder: (subst: (value: any) => string) => string) => {
    const { statement, vars } = prepare(statement_builder)
    return await this.pool.query(statement, vars)
  }
  
  listen = (cb: (evt: string, data: unknown) => void) => {
    const id = this.id++
    this.listeners[id] = cb
    return () => {delete this.listeners[id]}
  }
  protected notify = async (evt: string, data: unknown) => {
    for (const listener of Object.values(this.listeners)) {
      listener(evt, data)
    }
  }

  send = async (queue: string, work: unknown) => {
    await this.boss.send(queue, work as any)
  }

  work = async (queue: string, opts: unknown, cb: (work: unknown) => Promise<void>) => {
    await this.boss.work(queue, opts as any, cb)
    return () => {
      this.boss.stop().catch(error => console.error(error))
    }
  }
}

export class PgTable<T extends { id: Codec<string, string> }> implements DbTable<T> {
  constructor(public table: TypedSchema<T>, private db: PgDatabase) {}

  create = async (create: Create<T>)=> {
    const columns = Object.keys(this.table.field_codecs).filter(col => col in create.data) as Array<keyof T>
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
    const columns = Object.keys(this.table.field_codecs).filter(col => col in find.where) as Array<keyof T>
    if (columns.length <= 0) throw new Error('Missing where')
    const results = await this.db.raw(subst => `
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
    const results = await this.db.raw(subst => `
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
    const results = await this.db.raw(subst => `
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
    const results = await this.db.raw(subst => `
      with rows as (
        update ${JSON.stringify(this.table.name)}
        set ${set_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode((update.data as any)[key]))}`).join(', ')}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(update.where[key] as Decoded<T[keyof T]>))}`).join(' and ')};
        returning 1
      ) select count(*) as count from rows;
    `)
    return results.rows[0] as { count: number }
  }
  upsert = async (upsert: Upsert<T>) => {
    const insert_columns = Object.keys(this.table.field_codecs).filter(col => col in upsert.create) as Array<keyof T>
    if (insert_columns.length <= 0) throw new Error('Missing values')
    const update_columns = (upsert.update ? Object.keys(this.table.field_codecs).filter(col => upsert.update && col in upsert.update) : []) as Array<keyof T>
    dict.items(this.table.field_pk)
      .filter(({ value }) => value)
      .forEach(({ key }) => {
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
      on conflict (${dict.items(this.table.field_pk).filter(({ value }) => value).map(({ key }) => JSON.stringify(key)).join(',')})
      do update set ${update_columns.map(col => `${JSON.stringify(col)} = EXCLUDED.${JSON.stringify(col)}`).join(',')}
      returning *;
    `)
    if (results.rowCount !== 1) throw new Error('Expected one got none')
    return this.table.codec.decode(results.rows[0])
  }
  delete = async (delete_: Delete<T>) => {
    const where_columns = Object.keys(this.table.field_codecs).filter(col => col in delete_.where) as Array<keyof T>
    if (where_columns.length <= 0) throw new Error('Missing values')
    const results = await this.db.raw(subst => `
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
    const results = await this.db.raw(subst => `
      with rows as (
        delete from ${JSON.stringify(this.table.name)}
        where ${where_columns.map(key => `${JSON.stringify(key)} = ${subst(this.table.field_codecs[key].encode(delete_.where[key] as Decoded<T[keyof T]>))}`).join(' and ')};
        returning 1
      ) select count(*) as count from rows;
    `)
    return results.rows[0] as { count: number }
  }
}
