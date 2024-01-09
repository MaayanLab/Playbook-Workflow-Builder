import fetcher from "@/utils/next-rest-fetcher"
import type { FetchFileResponse } from "."

export async function clientFetchFile(file: { url: string }, session_id?: string): Promise<FetchFileResponse> {
  return await fetcher(`${session_id ? `/api/socket/${session_id}` : ''}/api/v1/components/core/file/fetch`, { headers: { 'Content-Type': 'application/json' }, method: 'POST', body: JSON.stringify(file) })
}
