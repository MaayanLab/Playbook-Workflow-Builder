import React from 'react'
import { fetcherPOST } from '@/utils/next-rest-fetcher'
import useAsyncEffect from 'use-async-effect'
import useSWRMutation from 'swr/mutation'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'

/**
 * Client interface for chatGPT -- checks if chatGPT is available and reports reason why it might not be and
 *  facilitates submitting a query.
 */
export function useChatGPT({ session_id }: { session_id?: string }) {
  const [chatGPTAvailable, setChatGPTAvailable] = React.useState(false)
  const { trigger, isMutating: isAugmentingWithChatGPT, error } = useSWRMutation(`${session_id ? `/api/socket/${session_id}` : ''}/api/chatgpt`, fetcherPOST)
  const { data: userSession } = useSessionWithId()
  const augmentWithChatGPT = React.useCallback((description: string) => {
    if (!chatGPTAvailable) throw new Error('ChatGPT is unavailable')
    return trigger({ description })
  }, [trigger, chatGPTAvailable])
  useAsyncEffect(async (isMounted) => {
    if (userSession?.user === null) return
    const req = await fetch(`/api/chatgpt`, { method: 'HEAD' })
    if (!isMounted()) return
    if (req.status === 200) setChatGPTAvailable(() => true)
  }, [userSession?.user])
  return {
    chatGPTAvailable,
    augmentWithChatGPT,
    isAugmentingWithChatGPT,
    errorAugmentingWithChatGPT:
      chatGPTAvailable ? error :
        !userSession?.user ? 'Please login'
        : 'ChatGPT is not configured',
  }
}
