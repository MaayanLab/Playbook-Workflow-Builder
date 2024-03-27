import type { UploadExampleFileResponse } from "."

export async function clientLoadExample(): Promise<UploadExampleFileResponse> {
  const res = await fetch(`/api/v1/components/data/gene_count_matrix/example.tsv`, { method: 'POST' })
  return await res.json()
}
