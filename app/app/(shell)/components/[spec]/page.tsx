import Client from '@/app/app/_clients/component'
import ClientOnly from '@/app/fragments/client-only'
import { decodeParam } from '@/app/app/_lib/params'

export default async function Page(props: { params: Promise<{ spec: string }> }) {
  const { spec } = await props.params
  return <ClientOnly><Client params={{ spec: decodeParam(spec) }} /></ClientOnly>
}
