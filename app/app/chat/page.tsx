import type { Metadata } from 'next'
import Client from '@/app/app/_clients/chat'
import ClientPage from '@/app/app/_clients/client-page'

export const metadata: Metadata = { title: 'Text to Workflow' }

export default function Page() {
  return <ClientPage><Client /></ClientPage>
}
