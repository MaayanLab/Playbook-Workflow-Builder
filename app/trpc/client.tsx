'use client'

import React from 'react'
import { QueryClient, QueryClientProvider } from '@tanstack/react-query'
import { createTRPCClient, httpBatchLink, httpSubscriptionLink, splitLink, TRPCClientError, type TRPCLink } from '@trpc/client'
import { createTRPCContext } from '@trpc/tanstack-react-query'
import type { TRPC_ERROR_CODE_KEY } from '@trpc/server/rpc'
import type { AppRouter } from './router'
import { TRPC_HTTP_PATH } from './constants'

// named exports (not destructured) so that they're visible to the client reference transform
const trpcContext = createTRPCContext<AppRouter>()
export const TRPCProvider = trpcContext.TRPCProvider
export const useTRPC = trpcContext.useTRPC
export const useTRPCClient = trpcContext.useTRPCClient

declare module '@trpc/client' {
  interface OperationContext {
    /**
     * Route this operation to the CAVATICA worker for this session
     *  (through the public server's tunnel) instead of this server.
     */
    session_id?: string
  }
}

/**
 * A link which delegates to a link constructed per session_id
 */
function sessionLink(factory: (session_id: string | undefined) => TRPCLink<AppRouter>): TRPCLink<AppRouter> {
  return (runtime) => {
    const links = new Map<string | undefined, ReturnType<TRPCLink<AppRouter>>>()
    return (props) => {
      const session_id = props.op.context.session_id
      let link = links.get(session_id)
      if (link === undefined) {
        link = factory(session_id)(runtime)
        links.set(session_id, link)
      }
      return link(props)
    }
  }
}

/** In a CAVATICA session, requests go to its worker through the public server's tunnel */
function trpcUrl(session_id: string | undefined) {
  return `${session_id ? `/api/socket/${session_id}` : ''}${TRPC_HTTP_PATH}`
}

function createClient() {
  return createTRPCClient<AppRouter>({
    links: [
      splitLink({
        // subscriptions are streamed as server-sent events
        condition: (op) => op.type === 'subscription',
        true: sessionLink((session_id) => httpSubscriptionLink({ url: trpcUrl(session_id) })),
        false: sessionLink((session_id) => httpBatchLink({ url: trpcUrl(session_id) })),
      }),
    ],
  })
}

export function TRPCReactProvider(props: React.PropsWithChildren<{}>) {
  const [queryClient] = React.useState(() => new QueryClient({
    defaultOptions: {
      queries: {
        // don't retry requests which will fail the same way again
        retry: (failureCount, error) => failureCount < 3 && !isTRPCClientError(error, ['UNAUTHORIZED', 'NOT_FOUND', 'FORBIDDEN', 'PAYMENT_REQUIRED', 'BAD_REQUEST']),
      },
    },
  }))
  const [trpcClient] = React.useState(createClient)
  return (
    <QueryClientProvider client={queryClient}>
      <TRPCProvider trpcClient={trpcClient} queryClient={queryClient}>
        {props.children}
      </TRPCProvider>
    </QueryClientProvider>
  )
}

export function isTRPCClientError(error: unknown, code?: TRPC_ERROR_CODE_KEY | TRPC_ERROR_CODE_KEY[]): error is TRPCClientError<AppRouter> {
  if (!(error instanceof TRPCClientError)) return false
  if (code === undefined) return true
  const codes = Array.isArray(code) ? code : [code]
  return codes.includes(error.data?.code)
}
