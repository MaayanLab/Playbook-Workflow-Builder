import fetcher from "@/utils/next-rest-fetcher"
import type { FetchFileResponse } from "."

export async function clientFetchFile(file: { url: string }): Promise<FetchFileResponse> {
  return await fetcher(`/api/v1/components/core/file/fetch`, { headers: { 'Content-Type': 'application/json' }, method: 'POST', body: JSON.stringify(file) })
}
