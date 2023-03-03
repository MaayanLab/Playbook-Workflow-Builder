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

export const config = {
  api: {
    bodyParser: false,
  },
}


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
  _path: z.array(z.string().transform(sanitize)).transform(ps => ps.join('/')),
}).transform(({_path}) => _path)

/**
 * Forward to component custom handler
 */
export default handler(async (req, res) => {
  const _path = QueryType.parse(req.query)
  let component_handlers: Record<string, RouteHandler> = {}
  if (!req.method) throw new UnsupportedMethodError()
  try {
    Object.assign(component_handlers, await require(`@/app/api/${_path}/route.ts`))
  } catch (e) {
    console.warn(e)
    throw new NotFoundError()
  }
  if (!(req.method in component_handlers)) throw new UnsupportedMethodError()
  await component_handlers[req.method](req, res, { throwOnError: false })
})
