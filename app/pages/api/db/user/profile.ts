import db from '@/app/db'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'

const BodyType = z.object({
  name: z.string().optional(),
  image: z.string().optional(),
  affiliation: z.string().optional(),
}).strict()

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) {
      res.status(401).end()
    } else if (req.method === 'GET') {
      const user = await db.objects.user.findUnique({ where: { id: session.user.id } })
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
    } else if (req.method === 'POST') {
      res.status(200).json(await db.objects.user.update({
        where: { id: session.user.id },
        data: BodyType.parse(JSON.parse(req.body)),
      }))
      
    } else {
      throw new Error('Unsupported method')
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
