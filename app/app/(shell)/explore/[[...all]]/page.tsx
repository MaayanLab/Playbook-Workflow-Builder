import type { Metadata } from 'next'
import Client from '@/app/app/_clients/explore'
import ClientOnly from '@/app/fragments/client-only'
import { decodeParam } from '@/app/app/_lib/params'

export const metadata: Metadata = { title: 'Playbook: Explorer' }

export default async function Page({ params }: { params: Promise<{ all?: string[] }> }) {
  const { all } = await params
  return <ClientOnly><Client initialPath={all?.map(decodeParam)} /></ClientOnly>
}
