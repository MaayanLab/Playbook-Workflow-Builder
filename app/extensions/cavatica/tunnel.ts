/**
 * The tunnel between this (public) server and CAVATICA workers.
 *
 * Workers can't be reached directly, instead they connect out to this server:
 *  - they subscribe to `tunnel.serve` (server-sent events) to receive requests
 *  - they respond with `tunnel.respond` mutations: a head, body chunks, then the end
 *
 * Requests to /api/socket/[session_id]/... are forwarded through the tunnel,
 *  responses are streamed back as they arrive (e.g. subscriptions).
 *
 * State is kept in memory: a session's worker & its requests must use the same
 *  ui instance (there is one).
 */
import { EventEmitter } from 'events'
import { randomUUID } from 'crypto'
import cache from '@/utils/global_cache'

export type TunnelRequest = {
  method: string
  path: string
  headers: Record<string, string>
  /** base64 encoded */
  body?: string
}

/** Sent to the worker */
export type WorkerEvent =
  | { type: 'request', id: string } & TunnelRequest
  | { type: 'abort', id: string }
  | { type: 'close' }

/** Sent by the worker */
export type WorkerResponse =
  | { type: 'head', id: string, status: number, headers: Record<string, string> }
  | { type: 'chunk', id: string, data: string }
  | { type: 'end', id: string }
  | { type: 'error', id: string, message: string }

export type ResponseHandlers = {
  onHead: (status: number, headers: Record<string, string>) => void
  onChunk: (data: Buffer) => void
  onEnd: () => void
  onError: (error: Error) => void
}

const workers = cache('tunnel:workers', () => new Map<string, (event: WorkerEvent) => void>())
const requests = cache('tunnel:requests', () => new Map<string, { session_id: string } & ResponseHandlers>())

/** `join:${session_id}` when a worker connects */
export const sessionEvents = cache('tunnel:events', () => new EventEmitter())

/** Register the worker for a session, returning a function to unregister it */
export function connectWorker(session_id: string, send: (event: WorkerEvent) => void) {
  workers.set(session_id, send)
  sessionEvents.emit(`join:${session_id}`)
  return () => {
    if (workers.get(session_id) === send) workers.delete(session_id)
  }
}

export function isWorkerConnected(session_id: string) {
  return workers.has(session_id)
}

/** Ask the session's worker to shut down */
export function closeSession(session_id: string) {
  workers.get(session_id)?.({ type: 'close' })
}

/** Forward a request to the session's worker, the returned function aborts it */
export function forward(session_id: string, request: TunnelRequest, handlers: ResponseHandlers) {
  const send = workers.get(session_id)
  if (!send) throw new Error('Worker not connected')
  const id = randomUUID()
  const done = () => { requests.delete(id) }
  requests.set(id, {
    session_id,
    onHead: handlers.onHead,
    onChunk: handlers.onChunk,
    onEnd: () => { done(); handlers.onEnd() },
    onError: (error) => { done(); handlers.onError(error) },
  })
  send({ type: 'request', id, ...request })
  return () => {
    if (requests.has(id)) {
      done()
      workers.get(session_id)?.({ type: 'abort', id })
    }
  }
}

/** Handle a response message from the worker of `session_id` */
export function respond(session_id: string, response: WorkerResponse) {
  const request = requests.get(response.id)
  // ignore responses to unknown (e.g. aborted) requests or another session's requests
  if (!request || request.session_id !== session_id) return
  if (response.type === 'head') request.onHead(response.status, response.headers)
  else if (response.type === 'chunk') request.onChunk(Buffer.from(response.data, 'base64'))
  else if (response.type === 'end') request.onEnd()
  else if (response.type === 'error') request.onError(new Error(response.message))
}
