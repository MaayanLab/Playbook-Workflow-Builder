/**
 * Expose component APIs
 */
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'

const QueryType = z.object({
  _path: z.array(z.string()),
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
 * Build a data structure for prefix-tree matching
 */
import { components } from '@/components'
import { create_prefix_tree_from_paths, search_prefix_tree } from '@/utils/prefix_tree'
const component_tree = create_prefix_tree_from_paths(components)

/**
 * Forward to component custom handler
 */
export default function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const { _path } = QueryType.parse(req.query)
    try {
      const { prefix: component, path } = search_prefix_tree(component_tree, _path.map(sanitize).join('/'))
      if (component === undefined || path === undefined) {
        throw new Error('Component not found')
      }
      const { default: handler } = require(`@/components/${component}/api/${path}`)
      return handler(req, res)
    } catch (e) {
      console.error(e)
      res.status(404).end()
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
