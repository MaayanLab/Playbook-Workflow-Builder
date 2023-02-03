import db from '@/app/db'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'

const QueryType = z.object({
  id: z.string(),
})

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (req.method !== 'GET') throw new Error('Unsupported method')
    const { id } = QueryType.parse(req.query)
    const user = await db.objects.user.findUnique({ where: { id } })
    if (user === null) {
      res.status(404).end()
    } else {
      res.status(200).json({
        image: user.image,
        name: user.name,
        email: user.email,
        affiliation: user.affiliation,
      })
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
