import type { Metadata } from 'next'
import Client from '@/app/app/_clients/community-playbooks'
import ClientPage from '@/app/app/_clients/client-page'

export const metadata: Metadata = { title: 'Community Playbooks' }

export default function Page() {
  return <ClientPage><Client /></ClientPage>
}
