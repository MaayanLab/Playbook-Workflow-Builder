import fs from 'fs'
import os from 'os'
import path from 'path'
import db from "@/app/db"
import { createHash } from 'crypto'
import type { Writable, Readable } from 'stream'
import type { SessionWithId } from "@/app/pages/api/auth/[...nextauth]"
import { toReadable } from '@/utils/readable'
import { fileAsStream } from '@/components/core/file/api/download'
import python from '@/utils/python'

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

export async function fileFromStream(reader_: Readable | ReadableStream<Uint8Array>, filename: string) {
  const reader = await toReadable(reader_)
  const tmp = path.join(fs.mkdtempSync(path.join(os.tmpdir(), 'ppwb-')), filename)
  const writer = fs.createWriteStream(tmp)
  const { sha256, size } = await statsFromStream(reader, writer)
  return { url: `file://${tmp.split(path.sep).join(path.posix.sep)}`, sha256, size, filename }
}

export type UploadFileResponse = Awaited<ReturnType<typeof uploadFile>>

export async function uploadFile(file: { url: string, size?: number, sha256?: string, filename: string }, session?: SessionWithId) {
  if (!file.sha256 || !file.size) {
    const stats = await statsFromStream(await fileAsStream(file))
    file.sha256 = stats.sha256
    if (file.size !== undefined && file.size !== stats.size) {
      console.warn('mismatched filesize')
    }
    file.size = stats.size
  }
  if (file.url.startsWith('file://')) {
    const origFile = {...file}
    file.url = `storage://${file.sha256}`
    await python('components.core.file.file_move', {
      kargs: [
        origFile,
        file
      ],
    })
  }
  const upload = await db.objects.upload.upsert({
    where: {
      url: file.url,
      sha256: file.sha256,
    },
    create: {
      url: file.url,
      sha256: file.sha256,
      size: file.size,
    },
  })
  const user_upload = await db.objects.user_upload.upsert({
    where: {
      user: session?.user?.id ?? null,
      upload: upload.id,
    },
    create: {
      user: session?.user?.id ?? null,
      upload: upload.id,
      filename: file.filename,
    },
  })
  return { url: `${(process.env.PUBLIC_URL||'').replace(/^https?:/, 'drs:')}/${user_upload.id}`, filename: file.filename, sha256: file.sha256, size: file.size }
}
