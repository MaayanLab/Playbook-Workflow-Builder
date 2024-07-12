import type { UploadExampleFileResponse } from "."

export async function clientLoadExample(): Promise<UploadExampleFileResponse> {
  const res = await fetch(`/api/v1/components/data/metadata_matrix/meta.tsv`, { method: 'POST' })
  return await res.json()
}
