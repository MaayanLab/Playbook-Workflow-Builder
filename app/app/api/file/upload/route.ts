import fs from 'fs/promises'
import os from 'os'
import path from 'path'
import { getRequestSessionWithId } from '@/app/extensions/auth/server'
import { uploadFile } from '@/components/core/file/server/upload'
import { routeResponse } from '@/app/app/api/_lib/response'
import { UnauthorizedError } from '@/spec/error'

export const runtime = 'nodejs'

export function POST(request: Request) {
  return routeResponse(async () => {
    const session = await getRequestSessionWithId(request)
    if (!session?.user) throw new UnauthorizedError()
    const form = await request.formData()
    const tmpDir = await fs.mkdtemp(path.join(os.tmpdir(), 'pwb-upload-'))
    try {
      const uploads: Record<string, unknown[]> = {}
      for (const [key, value] of form.entries()) {
        if (typeof value === 'string') continue
        const filename = value.name || 'upload'
        const filePath = path.join(tmpDir, `${crypto.randomUUID()}-${path.basename(filename)}`)
        await fs.writeFile(filePath, new Uint8Array(await value.arrayBuffer()))
        const uploaded = await uploadFile({
          url: `file://${filePath.split(path.sep).join(path.posix.sep)}`,
          size: value.size,
          filename,
        }, session)
        ;(uploads[key] ??= []).push(uploaded)
      }
      return Response.json(uploads)
    } finally {
      await fs.rm(tmpDir, { recursive: true, force: true })
    }
  })
}
