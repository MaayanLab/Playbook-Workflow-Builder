import { getServerSession } from 'next-auth/next'
import { authOptions } from '@/app/pages/api/auth/[...nextauth]'
import type { NextApiRequest, NextApiResponse } from 'next'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'

export async function getServerSessionWithId(req: NextApiRequest, res: NextApiResponse) {
  return await getServerSession(req, res, authOptions) as SessionWithId | null
}
