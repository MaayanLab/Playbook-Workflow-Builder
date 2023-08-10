import db from '@/app/db'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import handler from "@/utils/next-rest"

export default handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  if (req.method === 'GET') {
    const uploads = await db.objects.user_upload_complete.findMany({ where: { user: session.user.id } })
    uploads.forEach(upload => {
      upload.url = `${process.env.NEXT_PUBLIC_URL?.replace(/https?:\/\//, 'drs://')}/${upload.id}`
    })
    return res.status(200).json(uploads)
  } else {
    throw new UnsupportedMethodError()
  }
})
