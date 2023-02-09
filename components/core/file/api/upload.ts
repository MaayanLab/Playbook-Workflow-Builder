import multiparty from "multiparty"
import handler from '@/utils/next-rest'
import * as dict from '@/utils/dict'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError, UnsupportedMethodError } from "@/spec/error"
import db from "@/app/db"
import fs from 'fs'
import { createHash } from 'crypto'

function sha256FromFile(path: string) {
  return new Promise<string>((resolve, reject) => {
    const hash = createHash('sha256')
    fs.createReadStream(path)
      .on('error', (err) => reject(err))
      .on('data', (buf) => hash.update(buf))
      .on('end', () => resolve(hash.digest('hex')))
  })
}

export default handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  if (req.method !== 'POST') throw new UnsupportedMethodError()
  const form = new multiparty.Form({ autoFiles: true })
  const raw = await new Promise<{ fields: Record<string, string[]>, files: Record<string, multiparty.File[]> }>((resolve, reject) => {
    form.parse(req, function (err, fields, files) {
      if (err) reject({ err })
      else resolve({ fields, files })
    })
  })
  // TODO: store it in fsspec
  const arg_paths = dict.init(
    await Promise.all(
      dict.items(raw.files).map(async ({ key, value: files }) => {
        return {
          key,
          value: await Promise.all(files.map(async (file) => {
            const sha256 = await sha256FromFile(file.path)
            const upload = await db.objects.upload.upsert({
              where: {
                url: `file://${file.path}`,
                sha256,
              },
              create: {
                url: `file://${file.path}`,
                sha256,
                size: file.size,
              },
            })
            await db.objects.user_upload.upsert({
              where: {
                user: session.user.id,
                upload: upload.id,
              },
              create: {
                user: session.user.id,
                upload: upload.id,
                filename: file.originalFilename,
              },
            })
            return file.path
          })),
        }
      })
    )
  )
  res.status(200).json(arg_paths)
})
