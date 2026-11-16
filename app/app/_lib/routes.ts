/**
 * The routes of the playbook pages, parsed the same way by the server (pages validate them, e.g. to
 *  respond with a 404) and the browser (the views live in layouts, so that they persist while
 *  navigating within them, and read the route from the url).
 *
 * - /graph[/<graph_id>][/node/<node_id>][/extend|/suggest]
 * - /report[/<id>]
 * - either of those within /session/<session_id>
 */
import { decodeParam } from '@/app/app/_lib/params'

export type GraphRoute = {
  graph_id: string,
  /** the node(s) (comma separated), when given */
  node_id?: string,
  mode?: 'extend' | 'suggest',
}

export type ReportRoute = { id?: string }

export type PlaybookRoute =
  | { kind: 'graph', session_id?: string } & GraphRoute
  | { kind: 'report', session_id?: string } & ReportRoute
  | { kind: 'other', session_id?: string }

/** The path within /graph (already decoded), undefined if it isn't a valid graph route */
export function parseGraphPath(path: string[] = []): GraphRoute | undefined {
  let graph_id = 'start'
  let node_id: string | undefined
  let mode: GraphRoute['mode']
  const remaining = [...path]
  if (remaining[0] !== 'extend' && remaining[0] !== 'suggest') graph_id = remaining.shift() ?? 'start'
  if (remaining[0] === 'node') {
    remaining.shift()
    node_id = remaining.shift()
  }
  if (remaining[0] === 'extend' || remaining[0] === 'suggest') mode = remaining.shift() as GraphRoute['mode']
  if (remaining.length) return undefined
  return { graph_id, node_id, mode }
}

/** The path within /report (already decoded), undefined if it isn't a valid report route */
export function parseReportPath(path: string[] = []): ReportRoute | undefined {
  if (path.length > 1) return undefined
  return { id: path[0] }
}

/** The playbook route of a pathname (e.g. from usePathname, whose segments are still encoded) */
export function parsePathname(pathname: string): PlaybookRoute {
  const segments = pathname.split('/').filter(segment => segment !== '').map(decodeParam)
  let session_id: string | undefined
  if (segments[0] === 'session' && segments.length >= 2) {
    session_id = segments[1]
    segments.splice(0, 2)
  }
  if (segments[0] === 'graph') {
    const route = parseGraphPath(segments.slice(1))
    if (route) return { kind: 'graph', session_id, ...route }
  } else if (segments[0] === 'report') {
    const route = parseReportPath(segments.slice(1))
    if (route) return { kind: 'report', session_id, ...route }
  }
  return { kind: 'other', session_id }
}
