import useSWR from 'swr'

export default function UserIdentity(props: { user: string }) {
  const { data: user } = useSWR<{ image: string, name: string, email: string, affiliation: string }>(`/api/db/user/${props.user}/profile`)
  if (!user) return <span>a playbook partnership user</span>
  else return <span>{user.name} &lt;{user.email}&gt;{user.affiliation ? ` (${user.affiliation})` : ''}</span>
}
