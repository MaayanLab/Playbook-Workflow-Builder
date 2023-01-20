import Link from 'next/link'
import * as Auth from 'next-auth/react'

export default function Header() {
  const { data: session } = Auth.useSession()
  return (
    <div className="bg-primary flex flex-row justify-between place-items-center px-2">
      <Link href={process.env.NEXT_PUBLIC_LANDING_PAGE || '/'}>
        <h1 className="text-4xl font-bold p-2 cursor-pointer">Playbook Partnership Interactive Workflow Builder</h1>
      </Link>
      {session && session.user ?
        <span className="whitespace-nowrap">{session.user.email} (<button className="hover:underline" onClick={() => Auth.signOut()}>Sign Out</button>)</span>
        : <span className="whitespace-nowrap"><button className="hover:underline" onClick={() => Auth.signIn()}>Sign In</button></span>}
    </div>
  )
}
