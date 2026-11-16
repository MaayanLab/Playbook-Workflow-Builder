import type { Metadata } from 'next'
import Client from '@/app/app/_clients/workshop-2024-06-25'
import ClientPage from '@/app/app/_clients/client-page'

export const metadata: Metadata = { title: 'Playbook Workshop: June 25, 2024' }

export default function Page() {
  return <ClientPage><Client /></ClientPage>
}
