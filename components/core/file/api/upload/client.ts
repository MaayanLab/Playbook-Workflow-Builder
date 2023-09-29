import type { UploadFileResponse } from "."

export async function clientUploadFile(formData: FormData): Promise<Record<string, UploadFileResponse[]>> {
  const res = await fetch(`/api/v1/components/core/file/upload`, { method: 'POST', body: formData });
  return await res.json();
}
