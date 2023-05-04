import multiparty from "multiparty"
import handler from '@/utils/next-rest'
import * as dict from '@/utils/dict'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError } from "@/spec/error"
import { uploadFile } from "."

export const POST = handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user) throw new UnauthorizedError()
  const form = new multiparty.Form({ autoFiles: true })
  const raw = await new Promise<{ fields: Record<string, string[]>, files: Record<string, multiparty.File[]> }>((resolve, reject) => {
    form.parse(req, function (err, fields, files) {
      if (err) reject({ err })
      else resolve({ fields, files })
    })
  })
  const arg_paths = dict.init(
    await Promise.all(
      dict.items(raw.files).map(async ({ key, value: files }) => {
        return {
          key,
          value: await Promise.all(files.map(async (file) => {
            return await uploadFile({
              url: `file://${file.path}`,
              size: file.size,
              filename: file.originalFilename,
            })
          })),
        }
      })
    )
  )
  res.status(200).json(arg_paths)
})
