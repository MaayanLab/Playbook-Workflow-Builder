import React from 'react'
import { useMutation, useQuery } from '@tanstack/react-query'
import { useSessionWithId } from '@/app/extensions/auth/hooks'
import { useTRPC, useTRPCClient } from '@/app/trpc/client'

/**
 * Client interface for chatGPT -- checks if chatGPT is available and reports reason why it might not be and
 *  facilitates submitting a query.
 */
export function useChatGPT({ session_id }: { session_id?: string }) {
  const trpc = useTRPC()
  const trpcClient = useTRPCClient()
  const { data: userSession } = useSessionWithId()
  const { data: chatGPTAvailable = false } = useQuery({
    queryKey: ['chat-available', session_id],
    queryFn: () => trpcClient.chat.available.query(undefined, { context: { session_id } }),
    enabled: Boolean(userSession?.user),
    retry: false,
  })
  const augment = useMutation(trpc.chat.augment.mutationOptions({ trpc: { context: { session_id } } }))
  const augmentWithChatGPT = React.useCallback((description: string) => {
    if (!chatGPTAvailable) throw new Error('ChatGPT is unavailable')
    return augment.mutateAsync({ description })
  }, [augment, chatGPTAvailable])
  return {
    chatGPTAvailable,
    augmentWithChatGPT,
    isAugmentingWithChatGPT: augment.isPending,
    errorAugmentingWithChatGPT:
      chatGPTAvailable ? augment.error :
        !userSession?.user ? 'Please login'
        : 'ChatGPT is not configured',
  }
}
