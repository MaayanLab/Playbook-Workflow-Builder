import type { Session } from 'next-auth'

export default function UserAvatar({ session }: { session: Session | null }) {
  if (typeof session?.user?.image === 'string') {
    return <img src={session.user.image} />
  } else if (session?.user?.name) {
    const name_split = session.user.name.split(' ')
    const first_name = name_split[0]
    const last_name = name_split[name_split.length-1]
    return <span className="text-xl">{first_name[0].toUpperCase()}{last_name[0].toUpperCase()}</span>
  } else if (session?.user?.email) {
    return <span className="text-xl">{session.user.email[0].toUpperCase()}</span>
  } else {
    return <span className="text-xl">U</span>
  }
}