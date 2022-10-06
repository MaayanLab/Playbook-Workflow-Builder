/**
 * Expose component APIs
 */
import type { NextApiRequest, NextApiResponse } from 'next'

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
  const { _path } = req.query as { _path: string[] }
  const [component, ...path] = _path.map(sanitize)
  try {
    const { default: handler } = require(`@/components/${component}/api/${path.join('/')}`)
    return handler(req, res)
  } catch (e) {
    console.error(e)
    res.status(404).end()
  }
}
