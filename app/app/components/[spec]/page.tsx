import Client from '@/app/app/_clients/component'
import ClientPage from '@/app/app/_clients/client-page'
import { decodeParam } from '@/app/app/_lib/params'

export default async function Page(props: { params: Promise<{ spec: string }> }) {
  const { spec } = await props.params
  return <ClientPage><Client params={{ spec: decodeParam(spec) }} /></ClientPage>
}
