import type { Metadata } from 'next'
import Client from '@/app/app/_clients/account'
import ClientPage from '@/app/app/_clients/client-page'

export const metadata: Metadata = { title: 'Playbook: Account' }

export default function Page() {
  return <ClientPage><Client /></ClientPage>
}
