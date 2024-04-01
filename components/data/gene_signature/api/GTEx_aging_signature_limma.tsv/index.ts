import path from 'path'
import { uploadFile } from '@/components/core/file/api/upload'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'

export const exampleFile = path.resolve(
  process.env.APP_ROOT as string,
  'components',
  'data',
  'gene_signature',
  'api',
  'GTEx_aging_signature_limma.tsv',
  'GTEx_aging_signature_limma.tsv'
)

export type UploadExampleFileResponse = Awaited<ReturnType<typeof uploadExampleFile>>

export async function uploadExampleFile(session?: SessionWithId) {
  return await uploadFile({
    url: `${process.env.PUBLIC_URL}/api/v1/components/data/gene_signature/GTEx_aging_signature_limma.tsv`,
    filename: 'GTEx_aging_signature_limma.tsv',
  }, session)
}
