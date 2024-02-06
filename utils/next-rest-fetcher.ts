import { ResponseCodedError, NotFoundError, UnauthorizedError, UnsupportedMethodError, TimeoutError, UnboundError } from "@/spec/error"

export default async function fetcher<T>(...args: Parameters<typeof fetch>): Promise<T> {
  const req = await fetch(...args)
  try {
    if (req.status >= 200 && req.status < 300) {
      return await req.json()
    } else {
      const message = await req.text()
      if (req.status === 405) throw new UnsupportedMethodError(message)
      else if (req.status === 404) throw new NotFoundError(message)
      else if (req.status === 401) throw new UnauthorizedError(message)
      else if (req.status === 504) throw new TimeoutError(message)
      else if (req.status === 422) throw new UnboundError(message)
      else throw new ResponseCodedError(req.status, message)
    }
  } catch (e) {
    throw e
  }
}

export async function fetcherGET<T>(input: NodeJS.fetch.RequestInfo | URL) {
  return await fetcher<T>(input, { method: 'GET' })
}

export async function fetcherPOST<A, R>(input: NodeJS.fetch.RequestInfo | URL, { arg: body }: { arg?: A }) {
  return await fetcher<R>(input, {
    headers: {
      'Content-Type': 'application/json',
    },
    method: 'POST',
    body: body !== undefined ? JSON.stringify(body) : undefined,
  })
}
