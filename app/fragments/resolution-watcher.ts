/**
 * Watches the resolution of many processes over a single `fpprg.resolvedMany`
 *  subscription (instead of one per process), as each subscription holds a
 *  connection open and browsers limit concurrent HTTP/1.1 connections per host.
 *
 * - finished results are kept so that later watchers get them immediately
 * - the subscription only covers unfinished processes and is restarted (batched)
 *   when a new process is watched or one is refreshed
 */
import type { TRPCClient } from '@trpc/client'
import type { AppRouter } from '@/app/trpc/router'
import type { ResolvedLifecycle } from '@/app/trpc/routers/fpprg'

type Listener = (value: ResolvedLifecycle) => void

const isFinal = (value?: ResolvedLifecycle) => value !== undefined && value.type !== 'resolving'

export class ResolutionWatcher {
  private listeners = new Map<string, Set<Listener>>()
  private latest = new Map<string, ResolvedLifecycle>()
  private subscription: { unsubscribe: () => void } | undefined
  /** processes covered by the current subscription */
  private subscribed = new Set<string>()
  private scheduled: ReturnType<typeof setTimeout> | undefined
  private version = 0

  constructor(private client: TRPCClient<AppRouter>, private context: { session_id?: string }) {}

  watch = (id: string, listener: Listener) => {
    let listeners = this.listeners.get(id)
    if (listeners === undefined) this.listeners.set(id, listeners = new Set())
    listeners.add(listener)
    const value = this.latest.get(id)
    if (value !== undefined) listener(value)
    if (!isFinal(value) && !this.subscribed.has(id)) this.schedule()
    return () => {
      listeners.delete(listener)
      if (listeners.size > 0) return
      this.listeners.delete(id)
      if (this.listeners.size === 0) this.stop()
    }
  }

  /** Current lifecycle value for a process. Stable until the next event. */
  get = (id: string) => this.latest.get(id)

  /** Monotonic snapshot used when several process resolutions are consumed together. */
  getVersion = () => this.version

  /** Subscribe once to a set of process ids. */
  watchMany = (ids: string[], listener: () => void) => {
    const unwatch = ids.map(id => this.watch(id, listener))
    return () => unwatch.forEach(stop => stop())
  }

  /** Request the resolution again (e.g. after its output was deleted) */
  refresh = (id: string) => {
    this.latest.delete(id)
    this.subscribed.delete(id)
    if (this.listeners.has(id)) this.schedule()
  }

  private schedule() {
    // batch the watchers registered while a page renders into one subscription
    if (this.scheduled === undefined) {
      this.scheduled = setTimeout(() => {
        this.scheduled = undefined
        this.restart()
      }, 10)
    }
  }

  private restart() {
    const ids = [...this.listeners.keys()].filter(id => !isFinal(this.latest.get(id))).sort()
    this.subscription?.unsubscribe()
    this.subscription = undefined
    this.subscribed = new Set(ids)
    if (ids.length === 0) return
    this.subscription = this.client.fpprg.resolvedMany.subscribe(ids, {
      context: this.context,
      onData: ({ id, ...value }) => this.set(id, value),
      onError: (error) => {
        for (const id of ids) {
          if (!isFinal(this.latest.get(id))) this.set(id, { type: 'error', error: error.message })
        }
      },
    })
  }

  private set(id: string, value: ResolvedLifecycle) {
    this.latest.set(id, value)
    this.version += 1
    this.listeners.get(id)?.forEach(listener => listener(value))
  }

  private stop() {
    if (this.scheduled !== undefined) clearTimeout(this.scheduled)
    this.scheduled = undefined
    this.subscription?.unsubscribe()
    this.subscription = undefined
    this.subscribed.clear()
  }
}
