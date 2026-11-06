/**
 * Event bus & job queue used by the process resolution engine.
 *
 * - With postgres, events prefixed with `distributed:` go through LISTEN/NOTIFY
 *    (so they reach every process), insert notifications come from the
 *    `notify_trigger` triggers, and jobs are managed by pg-boss.
 * - With sqlite, everything happens in-process.
 */
import type * as pg from 'pg'
import PgBoss from 'pg-boss'
import createSubscriber, { type Subscriber } from 'pg-listen'
import { z } from 'zod'

export type Job = {
  id: string
  name: string
  priority: number
  state: string
  createdon: string
  completedon: string | null
}

export interface EventBus {
  listen: <T>(evt: string, cb: (data: T) => void) => () => void
  notify: <T>(evt: string, data: T) => void
  /** Announce that a row was inserted into `table`, emitted to listeners of `insert:${table}` */
  inserted: (table: string, id: string) => void
  /** Submit work to a queue, each id is processed by only one worker */
  send: (queue: string, work: { id: string, priority?: number }) => Promise<void>
  work: (queue: string, opts: { teamSize?: number, teamConcurrency?: number }, cb: (work: { data: { id: string } }) => Promise<void>) => Promise<() => void>
  jobs: () => Promise<Job[]>
}

class LocalListeners {
  private listeners: Record<string, Map<number, (data: unknown) => void>> = {}
  private id = 0
  add = (evt: string, cb: (data: unknown) => void) => {
    const id = this.id++
    if (!(evt in this.listeners)) this.listeners[evt] = new Map()
    this.listeners[evt].set(id, cb)
    return {
      first: this.listeners[evt].size === 1,
      remove: () => {
        this.listeners[evt]?.delete(id)
        if (this.listeners[evt]?.size === 0) {
          delete this.listeners[evt]
          return true
        }
        return false
      },
    }
  }
  emit = (evt: string, data: unknown) => {
    for (const listener of [...(this.listeners[evt]?.values() ?? [])]) {
      listener(data)
    }
  }
}

export class LocalEventBus implements EventBus {
  private listeners = new LocalListeners()
  private active: Record<string, Promise<void>> = {}

  listen = <T>(evt: string, cb: (data: T) => void) => this.listeners.add(evt, cb as (data: unknown) => void).remove as () => void
  notify = <T>(evt: string, data: T) => { this.listeners.emit(evt, data) }
  inserted = (table: string, id: string) => { this.notify(`insert:${table}`, { id }) }
  send = async (queue: string, work: { id: string, priority?: number }) => {
    this.notify(`boss:${queue}`, work)
  }
  work = async (queue: string, _opts: unknown, cb: (work: { data: { id: string } }) => Promise<void>) => {
    return this.listen(`boss:${queue}`, async (data) => {
      const { id } = z.object({ id: z.string() }).parse(data)
      if (!(id in this.active)) {
        this.active[id] = cb({ data: { id } })
          .catch(error => console.error(error))
          .finally(() => { delete this.active[id] })
      }
      return await this.active[id]
    })
  }
  jobs = async () => [] as Job[]
}

const InsertPayloadC = z.object({
  timestamp: z.string(),
  operation: z.string(),
  schema: z.string(),
  table: z.string(),
  label: z.array(z.string()),
  id: z.string(),
})

export class PgEventBus implements EventBus {
  private listeners = new LocalListeners()
  private subscriber: Subscriber
  private boss: PgBoss
  private _ready: Promise<void> | undefined

  constructor(connectionString: string, private pool: pg.Pool) {
    this.subscriber = createSubscriber({ connectionString })
    this.subscriber.events.on('error', (error) => {
      console.error('Fatal database connection error:', error)
      process.exit(1)
    })
    this.boss = new PgBoss(connectionString)
    this.boss.on('error', (error) => { console.error('PgBoss error:', error) })
  }

  /**
   * Connections are established on first use so that short-lived
   *  scripts (e.g. migrations) don't start pg-boss & a subscriber.
   */
  private ready = () => {
    if (this._ready === undefined) {
      this._ready = (async () => {
        process.on('exit', () => { this.subscriber.close() })
        await this.boss.start()
        await this.subscriber.connect()
      })()
      this._ready.catch((error) => {
        console.error('Failed to initialize subscriber', error)
        process.exit(1)
      })
      // relay row insertion triggers (see notify_trigger in the baseline migration)
      this.listen('distributed:on_insert', (rawPayload) => {
        const { table, operation, id } = InsertPayloadC.parse(rawPayload)
        if (operation === 'INSERT') this.listeners.emit(`insert:${table}`, { id })
      })
    }
    return this._ready
  }

  listen = <T>(evt: string, cb: (data: T) => void) => {
    const { first, remove } = this.listeners.add(evt, cb as (data: unknown) => void)
    const dbEvt = evt.startsWith('distributed:') ? evt.slice('distributed:'.length) : undefined
    if (evt.startsWith('insert:')) this.ready()
    if (first && dbEvt) {
      this.subscriber.notifications.on(dbEvt, (rawPayload) => this.listeners.emit(evt, rawPayload))
      this.ready().then(() => this.subscriber.listenTo(dbEvt)).catch(error => console.error(error))
    }
    return () => {
      if (remove() && dbEvt) {
        this.subscriber.notifications.removeAllListeners(dbEvt)
        this.ready().then(() => this.subscriber.unlisten(dbEvt)).catch(error => console.error(error))
      }
    }
  }

  /**
   * `distributed:` events are sent through the database and received by
   *  every listening process (including this one).
   */
  notify = <T>(evt: string, data: T) => {
    if (evt.startsWith('distributed:')) {
      const dbEvt = evt.slice('distributed:'.length)
      this.ready().then(() => this.subscriber.notify(dbEvt, data)).catch(error => console.error(error))
    } else {
      this.listeners.emit(evt, data)
    }
  }

  // insertions are announced by database triggers
  inserted = () => {}

  send = async (queue: string, work: { id: string, priority?: number }) => {
    await this.ready()
    await this.boss.send(queue, { id: work.id }, { singletonKey: work.id, priority: work.priority ?? 0 })
  }

  work = async (queue: string, opts: { teamSize?: number, teamConcurrency?: number }, cb: (work: { data: { id: string } }) => Promise<void>) => {
    await this.ready()
    await this.boss.work(queue, opts, cb as any)
    return () => { this.boss.stop().catch(error => console.error(error)) }
  }

  jobs = async () => {
    const results = await this.pool.query(`
      select name, priority, data, state, createdon, completedon
      from pgboss.job
      where data is not null
    `)
    return z.array(z.object({
      name: z.string(),
      priority: z.number(),
      data: z.object({ id: z.string() }),
      state: z.string(),
      createdon: z.date(),
      completedon: z.date().nullable(),
    })).parse(results.rows).map(({ data, createdon, completedon, ...job }) => ({
      ...job,
      id: data.id,
      createdon: createdon.toISOString(),
      completedon: completedon?.toISOString() ?? null,
    }))
  }
}
