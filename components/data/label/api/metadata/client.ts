import useSWR from 'swr'
import type { MetadataFromFileParams, MetadataFromFileResponse } from "."

export async function clientMetadataFromFile(file: MetadataFromFileParams[0]): Promise<MetadataFromFileResponse> {
  const res = await fetch(`/api/v1/components/data/label/metadata`, { method: 'POST', body: JSON.stringify(file) })
  return await res.json()
}

export function useClientMetadataFromFile(file: MetadataFromFileParams[0]) {
  return useSWR(file['url'], key => clientMetadataFromFile(file))
}
