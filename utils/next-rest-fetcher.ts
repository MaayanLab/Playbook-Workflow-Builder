import { ResponseCodedError, NotFoundError, UnauthorizedError, UnsupportedMethodError, TimeoutError, UnboundError } from "@/spec/error"

export default async function fetcher<T>(...args: Parameters<typeof fetch>): Promise<T> {
  const req = await fetch(...args)
  try {
    if (req.status >= 200 && req.status < 300) {
      return await req.json()
    } else {
      if (req.status === 405) throw new UnsupportedMethodError()
      else if (req.status === 404) throw new NotFoundError()
      else if (req.status === 401) throw new UnauthorizedError()
      else if (req.status === 504) throw new TimeoutError()
      else if (req.status === 422) throw new UnboundError()
      else throw new ResponseCodedError(req.status, await req.json())
    }
  } catch (e) {
    throw e
  }
}