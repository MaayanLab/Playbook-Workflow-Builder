/**
 * Expose APIs from @/app/api/*
 *
 * The reason for this redirection is to
 *   1. dispatch GET/POST/... to different handlers
 *   2. allow co-locating server route and client code to access that route
 */
import { z } from 'zod'
import handler, { RouteHandler } from '@/utils/next-rest'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import { components } from '@/components'
import { create_prefix_tree_from_paths, search_prefix_tree } from '@/utils/prefix_tree'
import APIRouter from '@/core/api/router'
import * as dict from '@/utils/dict'
import * as APIs from '@/app/api/server'

const api_router = new APIRouter()
dict.values(APIs).forEach(api => {api_router.add(api)})

export const config = {
  api: {
    bodyParser: false,
  },
}

/**
 * Build a data structure for prefix-tree matching
 */
const component_tree = create_prefix_tree_from_paths(components)

/**
 * Ensure we don't end up with weird strings
 */
function sanitize(component: string): string {
  const component_sanitized = component.replace(/[^A-Za-z0-9\._-]/g, '')
  if (/^\.*$/.exec(component_sanitized) !== null) {
    throw new Error('Invalid path component')
  }
  return component_sanitized
}

/**
 * Get sanitized path from query
 */
const QueryType = z.object({
  _path: z.array(z.string().transform(sanitize)),
}).transform(({_path}) => _path)

/**
 * Forward to component custom handler
 */
export default handler(async (req, res) => {
  const _path = QueryType.parse(req.query)
  let component_handlers: Record<string, RouteHandler> = {}
  if (!req.method) throw new UnsupportedMethodError()
  try {
    if (/^v\d+$/.exec(_path[0]) !== null && _path[1] === 'components') {
      const { prefix: component, path } = search_prefix_tree(component_tree, _path.slice(2).join('/'))
      if (component === undefined || path === undefined) throw new NotFoundError()
      Object.assign(component_handlers, await require(`@/components/${component}/api/${path}/route.ts`))
    } else {
      Object.assign(component_handlers, { GET: api_router.route, POST: api_router.route })
    }
  } catch (e) {
    console.warn(e)
    throw new NotFoundError()
  }
  if (!(req.method in component_handlers)) throw new UnsupportedMethodError(req.method)
  await component_handlers[req.method](req, res)
})
