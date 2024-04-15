import path from 'path'
import { uploadFile } from '@/components/core/file/api/upload'
import type { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'

export const exampleFile = path.resolve(
  process.env.APP_ROOT as string,
  'components',
  'data',
  'gene_matrix_transpose',
  'api',
  'Aging_Perturbations_from_GEO_Brain_up.gmt',
  'Aging_Perturbations_from_GEO_Brain_up.gmt'
)

export type UploadExampleFileResponse = Awaited<ReturnType<typeof uploadExampleFile>>

export async function uploadExampleFile(session?: SessionWithId) {
  return await uploadFile({
    url: `${process.env.PUBLIC_URL}/api/v1/components/data/gene_matrix_transpose/Aging_Perturbations_from_GEO_Brain_up.gmt`,
    filename: 'Aging_Perturbations_from_GEO_Brain_up.gmt',
  }, session)
}
