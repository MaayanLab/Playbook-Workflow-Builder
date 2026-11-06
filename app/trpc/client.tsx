import React from 'react'
import { QueryClient, QueryClientProvider } from '@tanstack/react-query'
import { createTRPCClient, createWSClient, httpBatchLink, splitLink, wsLink, TRPCClientError, type TRPCLink } from '@trpc/client'
import { createTRPCContext } from '@trpc/tanstack-react-query'
import type { TRPC_ERROR_CODE_KEY } from '@trpc/server/rpc'
import type { AppRouter } from './router'
import { TRPC_HTTP_PATH, TRPC_WS_PATH } from './constants'

export const { TRPCProvider, useTRPC, useTRPCClient } = createTRPCContext<AppRouter>()

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

function wsUrl(session_id: string | undefined) {
  const protocol = window.location.protocol === 'https:' ? 'wss:' : 'ws:'
  const search = session_id ? `?${new URLSearchParams({ session_id })}` : ''
  return `${protocol}//${window.location.host}${TRPC_WS_PATH}${search}`
}

function createClient() {
  return createTRPCClient<AppRouter>({
    links: [
      splitLink({
        condition: (op) => op.type === 'subscription',
        true: sessionLink((session_id) => wsLink({
          client: createWSClient({
            url: () => wsUrl(session_id),
            // only connect while there are active subscriptions
            lazy: { enabled: true, closeMs: 5_000 },
          }),
        })),
        false: sessionLink((session_id) => httpBatchLink({
          url: `${session_id ? `/api/socket/${session_id}` : ''}${TRPC_HTTP_PATH}`,
        })),
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
