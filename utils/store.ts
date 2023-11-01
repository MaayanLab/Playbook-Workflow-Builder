/**
 * This concept is taken from svelte, it's a convenient way
 *  for defining updateable data sources.
 */

export type Subscriber<T> = (value: T) => void
export type Unsubscriber = () => void
export type Readable<T> = {
  subscribe(run: Subscriber<T>): Unsubscriber
}
export type Writable<T> = Readable<T> & {
  set(value: T): void
}
export type StartStopNotifier<T> = (
  set: (value: T) => void,
) => void | (() => void)

export const writable = <T>(value?: T, start?: StartStopNotifier<T>): Writable<T> => {
  const listeners = new Set<Subscriber<T>>()
  const ctx = { value, unsub: undefined as void | (() => void) }
  return {
    set: (setValue) => {
      ctx.value = setValue
      listeners.forEach(listener => listener(setValue))
    },
    subscribe: (run) => {
      if (listeners.size === 0) {
        listeners.add(run)
        if (start !== undefined) {
          ctx.unsub = start(setValue => {
            ctx.value = setValue
            listeners.forEach(listener => listener(setValue))
          })
        }
      } else {
        listeners.add(run)
      }
      if (ctx.value !== undefined) {
        run(ctx.value)
      }
      return () => {
        listeners.delete(run)
        if (listeners.size === 0 && ctx.unsub !== undefined) {
          ctx.unsub()
        }
      }
    },
  }
}

export const readonly = <T>(store: Readable<T>): Readable<T> => {
  const { subscribe } = store
  return { subscribe }
}

export const readable = <T>(value?: T, start?: StartStopNotifier<T>): Readable<T> => readonly(writable(value, start))

export const readableFromPromise = <T>(value: Promise<T> | (() => Promise<T>)) =>
  readable<{ data?: T, error?: unknown }>({}, (set) => {
    (value instanceof Promise ? value : value())
      .then(data => set({ data }))
      .catch(error => set({ error }))
  })

export type Stores =
	| Readable<unknown>
	| [Readable<unknown>, ...Array<Readable<unknown>>]
	| Array<Readable<unknown>>;

/** One or more values from `Readable` stores. */
export type StoresValues<T> = T extends Readable<infer U>
	? U
	: { [K in keyof T]: T[K] extends Readable<infer U> ? U : never };

export const derived = <S extends Stores, T>(stores: S, fn: (values: StoresValues<S>) => T, initial_value?: T | undefined): Readable<T> =>
  readable(initial_value, (set) => {
    const stores_: Array<Readable<unknown>> = Array.isArray(stores) ? stores : [stores]
    const ctx = { ready: false }
    const store_ready = stores_.map(() => false)
    const values: unknown[] = stores_.map(() => undefined)
    const unsubs = stores_.map((store, i) =>
      store.subscribe(value => {
        values[i] = value
        if (!store_ready[i]) {
          store_ready[i] = true
          ctx.ready = !store_ready.some(value => value === false)
        }
        if (ctx.ready) {
          if (Array.isArray(stores)) {
            set(fn(values as StoresValues<S>))
          } else {
            set(fn(values[0] as StoresValues<S>))
          }
        }
      })
    )
    if (stores_.length === 0) {
      set(fn([] as StoresValues<S>))
    }
    return () => {
      unsubs.forEach(unsub => unsub())
    }
  })
