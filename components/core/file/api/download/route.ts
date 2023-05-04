import handler from '@/utils/next-rest'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError } from "@/spec/error"
import { z } from 'zod'
import { fileAsStream } from "."

const BodyType = z.object({
  url: z.string(),
})

export const POST = handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  const file = BodyType.parse(await req.json())
  const fileStream = await fileAsStream(file)
  fileStream.pipe(res)
})
