import fs from 'fs'
import os from 'os'
import path from 'path'
import db from "@/app/db"
import { createHash } from 'crypto'
import type { SessionWithId } from "@/app/pages/api/auth/[...nextauth]"
import type { Writable, Readable } from 'stream'

function statsFromStream(reader: Readable, writer?: Writable) {
  return new Promise<{ sha256: string, size: number }>((resolve, reject) => {
    let size = 0
    const hash = createHash('sha256')
    reader
      .on('error', (err) => reject(err))
      .on('data', (buf) => {
        size += buf.length
        hash.update(buf)
        if (writer) writer.write(buf, (err) => { if (err) reject(err) })
      })
      .on('end', () => {
        if (writer) writer.end()
        resolve({ sha256: hash.digest('hex'), size })
      })
  })
}

export async function fileFromStream(reader: Readable, originalFilename: string) {
  const tmp = path.join(fs.mkdtempSync(path.join(os.tmpdir(), 'ppwb-')), originalFilename)
  const writer = fs.createWriteStream(tmp)
  const { sha256, size } = await statsFromStream(reader, writer)
  return { path: tmp, sha256, size, originalFilename }
}

export type UploadFileResponse = Awaited<ReturnType<typeof uploadFile>>

export async function uploadFile(file: { path: string, size: number, sha256?: string, originalFilename: string }, session?: SessionWithId) {
  if (!file.sha256) {
    const stats = await statsFromStream(fs.createReadStream(file.path))
    file.sha256 = stats.sha256
    if (file.size !== stats.size) {
      file.size = stats.size
      console.warn('mismatched filesize')
    }
  }
  const url = `file://${file.path}`
  // TODO: store it in fsspec
  const upload = await db.objects.upload.upsert({
    where: {
      url,
      sha256: file.sha256,
    },
    create: {
      url,
      sha256: file.sha256,
      size: file.size,
    },
  })
  if (session && session.user && (session.user.id !== '00000000-0000-0000-0000-000000000000' || process.env.NODE_ENV === 'development')) {
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
  }
  return { url, filename: file.originalFilename, sha256: file.sha256, size: file.size }
}
