import db from "@/app/db"
import python from "@/utils/python"
import type { SessionWithId } from "@/app/pages/api/auth/[...nextauth]"

export type FetchFileResponse = Awaited<ReturnType<typeof fetchFile>>

export async function fetchFile({ url }: { url: string }, session?: SessionWithId) {
  const file: {
    url: string,
    filename: string,
    size: number,
    sha256: string,
  } = await python(
    'components.core.file.file_from_url', {
      kargs: [url]
    }
  )
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
