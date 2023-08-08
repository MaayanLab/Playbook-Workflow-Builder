import fetcher from "@/utils/next-rest-fetcher"
import type { UploadFileResponse } from "."

export async function clientUploadFile(formData: FormData): Promise<Record<string, UploadFileResponse[]>> {
  return await fetcher(`/api/v1/components/core/file/upload`, { method: 'POST', body: formData })
}
