import type { UploadExampleFileResponse } from "."

export async function clientLoadExample(): Promise<UploadExampleFileResponse> {
  const res = await fetch(`/api/v1/components/data/anndata/example.h5ad`, { method: 'POST' })
  return await res.json()
}
