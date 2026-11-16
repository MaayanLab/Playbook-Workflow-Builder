import { useQuery } from '@tanstack/react-query'
import { useTRPC } from '@/app/trpc/client'

export default function UserIdentity(props: { user: string }) {
  const trpc = useTRPC()
  const { data: user } = useQuery(trpc.users.profile.queryOptions({ id: props.user }))
  if (!user) return <span>a playbook workflow builder user</span>
  else return <span>{user.name}{user.email ? <>&lt;{user.email}&gt;</> : null}{user.affiliation ? ` (${user.affiliation})` : ''}</span>
}
