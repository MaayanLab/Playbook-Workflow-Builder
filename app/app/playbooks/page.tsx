import type { Metadata } from 'next'
import Client from '@/app/app/_clients/playbooks'
import ClientPage from '@/app/app/_clients/client-page'

export const metadata: Metadata = { title: 'Playbooks' }

export default function Page() {
  return <ClientPage><Client /></ClientPage>
}
