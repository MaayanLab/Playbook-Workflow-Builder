import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'

export type APIRoute<Q extends {} = any, O = unknown, B = any> =
  { path: string, method: 'GET', parameters: z.ZodType<Q>, call: (input: { query: Q }, req: NextApiRequest, res: NextApiResponse) => Promise<O> }
  | { path: string, method: 'POST', parameters: z.ZodType<Q>, requestBody: B, call: (input: { query: Q, body: B }, req: NextApiRequest, res: NextApiResponse) => Promise<O> }

export function API(path: string) {
  return ({
    query: <Q extends {}>(parameters: z.ZodType<Q>) => ({
      call: <O>(call: (input: { query: Q }, req: NextApiRequest, res: NextApiResponse) => Promise<O>) => ({
        build: () => ({ path, method: 'GET', parameters, call }) as APIRoute<Q, O>,
      }),
      body: <B>(requestBody: z.ZodType<B>) => ({
        call: <O>(call: (input: { query: Q, body: B }, req: NextApiRequest, res: NextApiResponse) => Promise<O>) => ({
          build: () => ({ path, method: 'POST', parameters, requestBody, call }) as APIRoute<Q, O, B>
        })
      }),
    }),
  })
}

export function APIInterface<T extends APIRoute<Q, O, B>, Q extends {} = any, O = unknown, B = any>(path: string, method: string) {
  return { path, method, call: z.custom<T['call']>() }
}

export type APIRouteInterface<Q extends {} = any, O = unknown, B = any> = ReturnType<typeof APIInterface<APIRoute<Q, O, B>, Q, O, B>>
