import React from 'react'
import fetcher from '@/utils/next-rest-fetcher'
import useAsyncEffect from 'use-async-effect'
import useSWRMutation from 'swr/mutation'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'

async function chatGPT(endpoint: string, { arg: description }: { arg: string }): Promise<string> {
  return await fetcher(endpoint, { method: 'POST', body: JSON.stringify({ description }) })
}

/**
 * Client interface for chatGPT -- checks if chatGPT is available and reports reason why it might not be and
 *  facilitates submitting a query.
 */
export function useChatGPT() {
  const [chatGPTAvailable, setChatGPTAvailable] = React.useState(false)
  const { trigger, isMutating: isAugmentingWithChatGPT, error } = useSWRMutation('/api/chatgpt', chatGPT)
  const { data: session } = useSessionWithId()
  const augmentWithChatGPT = React.useCallback((description: string) => {
    if (!chatGPTAvailable) throw new Error('ChatGPT is unavailable')
    return trigger(description)
  }, [trigger, chatGPTAvailable])
  useAsyncEffect(async (isMounted) => {
    if (session?.user === null) return
    const req = await fetch('/api/chatgpt', { method: 'HEAD' })
    if (!isMounted()) return
    if (req.status === 200) setChatGPTAvailable(() => true)
  }, [session?.user])
  return {
    chatGPTAvailable,
    augmentWithChatGPT,
    isAugmentingWithChatGPT,
    errorAugmentingWithChatGPT:
      chatGPTAvailable ? error :
        !session?.user ? 'Please login'
        : 'ChatGPT is not configured',
  }
}
