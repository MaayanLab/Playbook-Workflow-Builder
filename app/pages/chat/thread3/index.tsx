import React from 'react'
import { GPTAssistantThread3Create } from "@/app/api/client"
import { useAPIMutation } from "@/core/api/client"
import { useRouter } from "next/router"

export default function ChatThread() {
  const router = useRouter()
  const {trigger} = useAPIMutation(GPTAssistantThread3Create)
  React.useEffect(() => {
    (async () => {
      const thread_id = await trigger()
      router.push(`/chat/thread3/${thread_id}`)
    })()
  }, [])
  return (
    <div className="alert alert-info">
      Starting...
    </div>
  )
}