/**
 * Expose component APIs
 */
import type { NextApiRequest, NextApiResponse } from 'next'
import * as t from 'io-ts'
import decodeOrThrow from '@/utils/decodeOrThrow'

const QueryType = t.type({
  _path: t.array(t.string),
})

export const config = {
  api: {
    bodyParser: false,
  },
}

/**
 * Ensure we don't end up with weird strings
 */
function sanitize(component: string): string {
  return component.replace(/[^A-Za-z0-9_-]/g, '')
}

/**
 * Forward to component custom handler
 */
export default function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const { _path } = decodeOrThrow(QueryType, req.query)
    const [component, ...path] = _path.map(sanitize)
    try {
      const { default: handler } = require(`@/components/${component}/api/${path.join('/')}`)
      return handler(req, res)
    } catch (e) {
      console.error(e)
      res.status(404).end()
    }
  } catch (e) {
    console.error(e)
    res.status(500).end(e.toString())
  }
}
