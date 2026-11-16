import type { Metadata } from 'next'
import Client from '@/app/app/_clients/sitemap'
import ClientOnly from '@/app/fragments/client-only'

export const metadata: Metadata = { title: 'Playbook: Site Map' }

export default function Page() {
  return <ClientOnly><Client /></ClientOnly>
}
