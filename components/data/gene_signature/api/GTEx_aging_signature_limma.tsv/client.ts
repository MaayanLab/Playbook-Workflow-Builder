import type { UploadExampleFileResponse } from "."

export async function clientLoadExample(): Promise<UploadExampleFileResponse> {
  const res = await fetch(`/api/v1/components/data/gene_signature/GTEx_aging_signature_limma.tsv`, { method: 'POST' })
  return await res.json()
}
