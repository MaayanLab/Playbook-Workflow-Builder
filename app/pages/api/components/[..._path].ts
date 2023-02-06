/**
 * Expose component APIs
 */
import { z } from 'zod'
import handler from '@/utils/next-rest'
import { NotFoundError } from '@/spec/error'
import type { NextApiRequest, NextApiResponse } from 'next'

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
export default handler(async (req, res) => {
  const { _path } = QueryType.parse(req.query)
  const { prefix: component, path } = search_prefix_tree(component_tree, _path.map(sanitize).join('/'))
  if (component === undefined || path === undefined) throw new NotFoundError()
  let handler: (req: NextApiRequest, res: NextApiResponse) => Promise<void>
  try {
    handler = require(`@/components/${component}/api/${path}`).default
  } catch (e) {
    console.error(e)
    throw new NotFoundError()
  }
  await handler(req, res)
})
