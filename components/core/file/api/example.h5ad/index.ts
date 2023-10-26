import path from 'path'
import { uploadFile } from '@/components/core/file/api/upload'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'

export const exampleFile = path.resolve(
  process.env.APP_ROOT as string,
  'components',
  'core',
  'file',
  'api',
  'example.h5ad',
  'example.h5ad'
)

export type UploadExampleFileResponse = Awaited<ReturnType<typeof uploadExampleFile>>

export async function uploadExampleFile(session?: SessionWithId) {
  return await uploadFile({
    url: `file://${exampleFile.split(path.sep).join(path.posix.sep)}`,
    filename: 'example.h5ad',
  }, session)
}
