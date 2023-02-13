import Link from 'next/link'
import * as Auth from 'next-auth/react'
import type { Session } from 'next-auth'
import { useRuntimeConfig } from '@/app/fragments/config'

function UserAvatar({ session }: { session: Session | null }) {
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

export default function Header() {
  const { data: session } = Auth.useSession()
  const runtimeConfig = useRuntimeConfig()
  return (
    <div className="bg-primary flex flex-row justify-between place-items-center px-2">
      <Link href={runtimeConfig.NEXT_PUBLIC_LANDING_PAGE}>
        <h1 className="text-4xl font-bold p-2 cursor-pointer">Playbook Partnership Interactive Workflow Builder</h1>
      </Link>
      {session && session.user ?
        <div className="dropdown dropdown-end">
          <label tabIndex={0} className="btn btn-ghost btn-circle avatar placeholder">
            <div className="bg-neutral-focus text-neutral-content rounded-full w-16">
              <UserAvatar session={session} />
            </div>
          </label>
          <ul tabIndex={0} className="menu menu-compact dropdown-content mt-3 p-2 shadow bg-base-100 rounded-box w-52">
            <li><Link href="/account">Settings</Link></li>
            <li><Link href="/api/auth/signout">Sign Out</Link></li>
          </ul>
        </div>
        : <Link href="/api/auth/signin"><button className="btn btn-ghost">Sign in</button></Link>}
    </div>
  )
}
