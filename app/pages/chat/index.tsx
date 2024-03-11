import React from 'react'
import { GPTAssistantCreate } from "@/app/api/client"
import { useAPIMutation } from "@/core/api/client"
import { useRouter } from "next/router"
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import { signIn } from 'next-auth/react'

export default function ChatThread() {
  const router = useRouter()
  const auth = useSessionWithId({ required: true })
  const {trigger} = useAPIMutation(GPTAssistantCreate)
  React.useEffect(() => {
    if (!auth.data?.user?.id) {
      signIn()
      return
    }
    (async () => {
      const thread_id = await trigger()
      router.push(`/chat/${thread_id}`)
    })()
  }, [auth])
  return (
    <div className="alert alert-info">
      Starting...
    </div>
  )
}