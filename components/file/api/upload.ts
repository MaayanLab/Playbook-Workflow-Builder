import multiparty from "multiparty"
import type { NextApiRequest, NextApiResponse } from 'next'

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  const form = new multiparty.Form({ autoFiles: true })
  const raw = await new Promise<{ fields: Record<string, string[]>, files: Record<string, multiparty.File[]> }>((resolve, reject) => {
    form.parse(req, function (err, fields, files) {
      if (err) reject({ err })
      else resolve({ fields, files })
    })
  })
  // TODO: store it in fsspec
  const paths = Object.keys(raw.files)
    .reduce((paths, key) =>
      ({ ...paths, [key]: raw.files[key].map(file => file.path) }),
      {} as Record<string, string[]>
    )
  res.status(200).json(paths)
}
