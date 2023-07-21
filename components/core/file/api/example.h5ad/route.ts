import handler from '@/utils/next-rest'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError } from "@/spec/error"
import fs from 'fs'
import { exampleFile, uploadExampleFile } from '.'

export const GET = handler(async (req, res) => {
  fs.createReadStream(exampleFile).pipe(res)
})

export const POST = handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  res.status(200).json(await uploadExampleFile(session)
  )
})
