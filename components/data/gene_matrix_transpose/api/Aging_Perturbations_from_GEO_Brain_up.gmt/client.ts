import type { UploadExampleFileResponse } from "."

export async function clientLoadExample(): Promise<UploadExampleFileResponse> {
  const res = await fetch(`/api/v1/components/data/gene_matrix_transpose/Aging_Perturbations_from_GEO_Brain_up.gmt`, { method: 'POST' })
  return await res.json()
}
