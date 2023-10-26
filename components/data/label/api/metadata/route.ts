import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError } from "@/spec/error"
import handler from "@/utils/next-rest"
import getRawBody from 'raw-body'
import { metadataFromFile } from "."

export const POST = handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  const rawBody = (await getRawBody(req)).toString()
  const body = JSON.parse(rawBody)
  const metadata = await metadataFromFile(body)
  res.json(metadata)
})
