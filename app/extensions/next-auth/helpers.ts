import { getServerSession } from 'next-auth/next'
import { SYSTEM_USER_ID } from '@/db/constants'
import { authOptions } from '@/app/pages/api/auth/[...nextauth]'
import type { NextApiRequest, NextApiResponse } from 'next'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'

export async function getServerSessionWithId(req: NextApiRequest, res: NextApiResponse): Promise<SessionWithId | null> {
  if (req.headers.authorization === `Token ${process.env.NEXTAUTH_SECRET}`) {
    // worker authentication
    // TODO: workers should work on behalf of the user who triggered them
    return {
      user: {
        id: SYSTEM_USER_ID
      },
      expires: (new Date(Date.now() + 1000*60*60*24)).toISOString(),
    } as SessionWithId
  }
  return await getServerSession(req, res, authOptions)
}
