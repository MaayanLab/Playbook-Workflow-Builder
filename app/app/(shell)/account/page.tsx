import type { Metadata } from 'next'
import Client from '@/app/app/_clients/account'
import ClientOnly from '@/app/fragments/client-only'

export const metadata: Metadata = { title: 'Playbook: Account' }

export default function Page() {
  return <ClientOnly><Client /></ClientOnly>
}
