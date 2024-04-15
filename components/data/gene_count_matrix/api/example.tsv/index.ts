import path from 'path'
import { uploadFile } from '@/components/core/file/api/upload'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'

export const exampleFile = path.resolve(
  process.env.APP_ROOT as string,
  'components',
  'data',
  'gene_count_matrix',
  'api',
  'example.tsv',
  'example.tsv'
)

export type UploadExampleFileResponse = Awaited<ReturnType<typeof uploadExampleFile>>

export async function uploadExampleFile(session?: SessionWithId) {
  return await uploadFile({
    url: `${process.env.PUBLIC_URL}/api/v1/components/data/gene_count_matrix/example.tsv`,
    filename: 'example.tsv',
  }, session)
}
