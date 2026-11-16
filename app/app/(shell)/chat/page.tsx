import type { Metadata } from 'next'
import Client from '@/app/app/_clients/chat'
import ClientOnly from '@/app/fragments/client-only'

export const metadata: Metadata = { title: 'Text to Workflow' }

export default function Page() {
  return <ClientOnly><Client /></ClientOnly>
}
