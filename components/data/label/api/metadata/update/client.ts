import useSWRMutation from 'swr/mutation'
import type { UpdateMetadataParams, UpdateMetadataResponse } from "."

export async function clientUpdateMetadata(req: UpdateMetadataParams[0]): Promise<UpdateMetadataResponse> {
  const res = await fetch(`/api/v1/components/data/label/metadata/update`, { method: 'POST', body: JSON.stringify(req) })
  return await res.json()
}

export function useClientUpdateMetadata(key: string) {
  return useSWRMutation(key, (key, { arg }: { arg: UpdateMetadataParams[0] }) => clientUpdateMetadata(arg))
}
