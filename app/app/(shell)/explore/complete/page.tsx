import Client from '@/app/app/_clients/explore-complete'
import ClientOnly from '@/app/fragments/client-only'

export default function Page() {
  return <ClientOnly><Client /></ClientOnly>
}
