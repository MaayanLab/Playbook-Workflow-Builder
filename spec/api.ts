import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'

export type APIRoute<Q extends {} = any, O = unknown, B extends z.ZodType | FormData = any> = 
  { path: string, method: 'GET', parameters: z.ZodType<Q>, call: (input: { query: Q }, req: NextApiRequest, res: NextApiResponse) => Promise<O> }
  | { path: string, method: 'POST', parameters: z.ZodType<Q>, requestBody: B extends z.ZodType<infer T> ? T : B, call: (input: { query: Q, body: B extends z.ZodType<infer T> ? T : B }, req: NextApiRequest, res: NextApiResponse) => Promise<O> }

export type APIRouteInterface<Q extends {} = any, O = unknown, B extends z.ZodType | FormData = any> = 
  { path: string, method: 'GET', call: z.ZodType<(input: { query: Q }, req: NextApiRequest, res: NextApiResponse) => Promise<O>> }
  | { path: string, method: 'POST', call: z.ZodType<(input: { query: Q, body: B extends z.ZodType<infer T> ? T : B }, req: NextApiRequest, res: NextApiResponse) => Promise<O>> }

export function API(path: string) {
  return ({
    query: <Q extends {}>(parameters: z.ZodType<Q>) => ({
      call: <O>(call: (input: { query: Q }, req: NextApiRequest, res: NextApiResponse) => Promise<O>) => ({
        build: () => ({ path, method: 'GET', parameters, call }) as APIRoute<Q, O>,
      }),
      body: <B extends z.ZodType | FormData>(requestBody: B) => ({
        call: <O>(call: (input: { query: Q, body: B extends z.ZodType<infer T> ? T : B }, req: NextApiRequest, res: NextApiResponse) => Promise<O>) => ({
          build: () => ({ path, method: 'POST', parameters, requestBody, call }) as APIRoute<Q, O, B>
        })
      }),
    }),
  })
}
