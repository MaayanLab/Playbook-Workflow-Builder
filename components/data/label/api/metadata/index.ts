import type { SessionWithId } from "@/app/pages/api/auth/[...nextauth]"
import python from "@/utils/python"

export type MetadataFromFileParams = Parameters<typeof metadataFromFile>
export type MetadataFromFileResponse = Awaited<ReturnType<typeof metadataFromFile>>

export async function metadataFromFile(file: { url: string, size?: number, sha256?: string, filename: string }, session?: SessionWithId) {
  return await python('components.data.label.get_metadata_from_anndata', {
    kargs: [file],
  }) as Record<string, Record<string, string>>
}
