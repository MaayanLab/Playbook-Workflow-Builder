import type { UploadExampleFileResponse } from "."

export async function clientLoadExample(): Promise<UploadExampleFileResponse> {
  const res = await fetch(`/api/v1/components/gly_gen/glycosight/`, { method: 'POST' })
  return await res.json()
}
