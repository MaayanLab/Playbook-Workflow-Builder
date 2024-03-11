import React from 'react'
import { GPTAssistantThreadCreate } from "@/app/api/client"
import { useAPIMutation } from "@/core/api/client"
import { useRouter } from "next/router"

export default function ChatThread() {
  const router = useRouter()
  const {trigger} = useAPIMutation(GPTAssistantThreadCreate)
  React.useEffect(() => {
    (async () => {
      const thread_id = await trigger()
      router.push(`/chat/thread/${thread_id}`)
    })()
  }, [])
  return (
    <div className="alert alert-info">
      Starting...
    </div>
  )
}