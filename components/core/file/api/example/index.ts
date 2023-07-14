import path from 'path'
import { uploadFile } from '@/components/core/file/api/upload'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'

export const exampleFile = path.resolve(
  process.cwd(),
  '..',
  'components',
  'core',
  'file',
  'api',
  'example',
  'example_matrix.tsv'
)

export type UploadExampleFileResponse = Awaited<ReturnType<typeof uploadExampleFile>>

export async function uploadExampleFile(session?: SessionWithId) {
  return await uploadFile({
    url: `file://${exampleFile.split(path.sep).join(path.posix.sep)}`,
    filename: 'example_matrix.tsv',
  }, session)
}
