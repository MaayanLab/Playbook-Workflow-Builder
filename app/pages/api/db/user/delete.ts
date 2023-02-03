import db from '@/app/db'
import type { NextApiRequest, NextApiResponse } from 'next'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (req.method !== 'POST') throw new Error('Unsupported method')
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) {
      res.status(401).end()
    } else {
      const deletedUser = await db.objects.user.delete({ where: { id: session.user.id } })
      if (deletedUser === null) {
        res.status(404).end()
      } else {
        res.status(200).json(deletedUser)
      }
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
