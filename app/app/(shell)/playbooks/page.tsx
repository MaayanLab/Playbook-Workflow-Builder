import type { Metadata } from 'next'
import Client from '@/app/app/_clients/playbooks'
import ClientOnly from '@/app/fragments/client-only'

export const metadata: Metadata = { title: 'Playbooks' }

export default function Page() {
  return <ClientOnly><Client /></ClientOnly>
}
