/**
 * Promise helpers in the spirit of `Promise.all`, for when running everything at once is too much
 *  (e.g. memory): the tasks are given as functions so that they're only started when allowed.
 */

/**
 * A limit of `max` tasks running at once: `run(task)` starts `task` once fewer than `max` of the
 *  tasks given to `run` are running (in the order they're given), and resolves to its result.
 *  Share one `run` to limit tasks started in different places (e.g. across requests).
 */
export function limit(max: number) {
  let available = Math.max(1, max)
  const waiting: (() => void)[] = []
  return async function run<T>(task: () => Promise<T>): Promise<T> {
    if (available > 0) available--
    else await new Promise<void>(resolve => waiting.push(resolve))
    try {
      return await task()
    } finally {
      // hand the slot straight to the next task, if there is one
      const next = waiting.shift()
      if (next) next()
      else available++
    }
  }
}

/**
 * Like `Promise.all`, with at most `max` of the `tasks` running at once.
 *  The results are in the order of the tasks; it rejects with the first failure
 *  (the tasks already running finish, the others aren't started).
 */
export async function allConcurrent<T>(tasks: (() => Promise<T>)[], max: number): Promise<T[]> {
  const run = limit(max)
  let failed = false
  return await Promise.all(tasks.map(task => run(async () => {
    if (failed) throw new Error('not started: an earlier task failed')
    try {
      return await task()
    } catch (error) {
      failed = true
      throw error
    }
  })))
}

/**
 * Like `Promise.all`, running the `tasks` one at a time, in order.
 */
export function allSequential<T>(tasks: (() => Promise<T>)[]): Promise<T[]> {
  return allConcurrent(tasks, 1)
}
