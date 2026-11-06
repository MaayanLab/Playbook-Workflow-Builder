import fetcher from "@/utils/next-rest-fetcher"
import type { UploadFileResponse } from "../server/upload"

export async function clientUploadFile(formData: FormData, session_id?: string): Promise<Record<string, UploadFileResponse[]>> {
  return await fetcher(`${session_id ? `/api/socket/${session_id}` : ''}/api/file/upload`, { method: 'POST', body: formData })
}
