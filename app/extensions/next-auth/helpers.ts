import { getServerSession } from 'next-auth/next'
import { authOptions } from '@/app/pages/api/auth/[...nextauth]'
import type { NextApiRequest, NextApiResponse } from 'next'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'

export async function getServerSessionWithId(req: NextApiRequest, res: NextApiResponse) {
  if (req.headers.authorization === `Token ${process.env.NEXTAUTH_SECRET}`) {
    // worker authentication
    // TODO: workers should work on behalf of the user who triggered them
    return {
      user: {
        id: '00000000-0000-0000-0000-000000000000'
      }
    }
  }
  return await getServerSession(req, res, authOptions) as SessionWithId | null
}
