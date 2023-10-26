import { z } from 'zod'
import handler from '@/utils/next-rest'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError } from "@/spec/error"
import { fetchFile } from "."
import getRawBody from 'raw-body'

const BodyType = z.object({
  url: z.string(),
})

export const POST = handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  const body = (await getRawBody(req)).toString()
  const file = BodyType.parse(JSON.parse(body))
  const result = await fetchFile(file, session)
  res.status(200).json(result)
})
