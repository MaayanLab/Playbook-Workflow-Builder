/**
 * `fetch`, retrying requests which are rate limited (429 Too Many Requests) once the server says
 *  they may be retried: e.g. Enrichr allows one request per second, and the workers can run steps
 *  which use it at the same time.
 */

const sleep = (ms: number, signal?: AbortSignal | null) => new Promise<void>((resolve, reject) => {
  if (signal?.aborted) return reject(signal.reason)
  const timeout = setTimeout(() => { signal?.removeEventListener('abort', onAbort); resolve() }, ms)
  const onAbort = () => { clearTimeout(timeout); reject(signal?.reason) }
  signal?.addEventListener('abort', onAbort, { once: true })
})

/**
 * How long the response says to wait before retrying, in ms: `X-Retry-In` (e.g. `996.99ms`, which
 *  Enrichr sends, more precise than) `Retry-After` (seconds, or an HTTP date).
 */
export function retryDelay(response: Response): number | undefined {
  const retryIn = /^\s*([\d.]+)\s*(ms|s)?\s*$/.exec(response.headers.get('x-retry-in') ?? '')
  if (retryIn) return +retryIn[1] * (retryIn[2] === 's' ? 1000 : 1)
  const retryAfter = response.headers.get('retry-after')?.trim()
  if (!retryAfter) return undefined
  if (/^\d+(\.\d+)?$/.test(retryAfter)) return +retryAfter * 1000
  const date = Date.parse(retryAfter)
  return Number.isNaN(date) ? undefined : Math.max(0, date - Date.now())
}

export type RetryOptions = {
  /** how many times a rate limited request is retried (default 10) */
  retries?: number,
  /** the longest the server may have us wait (default 60s), beyond which the 429 is returned */
  maxDelay?: number,
}

/**
 * Like `fetch`, but a 429 response is retried after the delay it gives (or with exponential
 *  backoff if it gives none), up to `retries` times; any other response is returned as it is.
 *  (the request's body must be one which can be sent again, e.g. a string or FormData, not a stream)
 */
export async function fetchRetrying(input: string | URL, init?: RequestInit, { retries = 10, maxDelay = 60_000 }: RetryOptions = {}) {
  for (let attempt = 0; ; attempt++) {
    const response = await fetch(input, init)
    if (response.status !== 429 || attempt >= retries) return response
    const delay = retryDelay(response) ?? 1000 * 2 ** attempt
    if (delay > maxDelay) return response
    // (the body isn't needed, release the connection)
    await response.body?.cancel().catch(() => {})
    // with some jitter, so requests which were limited together don't retry together
    await sleep(delay + Math.random() * 250 * (attempt + 1), init?.signal)
  }
}
