import type { SessionWithId } from "@/app/pages/api/auth/[...nextauth]"
import python from "@/utils/python"

export type UpdateMetadataParams = Parameters<typeof updateMetadataColumn>
export type UpdateMetadataResponse = Awaited<ReturnType<typeof updateMetadataColumn>>

export async function updateMetadataColumn(props: {
  file: { url: string, size?: number, sha256?: string, filename: string, description?: string }
  data: Record<string, Record<string, string>>,
}, session?: SessionWithId) {
  return await python('components.data.label.update_anndata_metadata', {
    kargs: [props.file, props.data],
  }) as any
}
