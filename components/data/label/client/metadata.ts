import { useQuery } from '@tanstack/react-query'
import type { MetadataFromFileParams } from "../server/metadata"
import { useTRPC } from '@/app/trpc/client'

export function useClientMetadataFromFile(file: MetadataFromFileParams[0]) {
  const trpc = useTRPC()
  return useQuery(trpc.components.label.metadata.queryOptions(file))
}
