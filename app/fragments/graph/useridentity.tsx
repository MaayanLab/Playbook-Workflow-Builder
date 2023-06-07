import useSWR from 'swr'
import { fetcherGET } from '@/utils/next-rest-fetcher'

export default function UserIdentity(props: { user: string }) {
  const { data: user } = useSWR<{ image: string, name: string, email: string, affiliation: string }>(() => `/api/db/user/${props.user}/profile`, fetcherGET)
  if (!user) return <span>a playbook partnership user</span>
  else return <span>{user.name}{user.email ? <>&lt;{user.email}&gt;</> : null}{user.affiliation ? ` (${user.affiliation})` : ''}</span>
}
