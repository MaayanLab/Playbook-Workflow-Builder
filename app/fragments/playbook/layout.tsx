import Link from 'next/link'
import * as Auth from 'next-auth/react'
import type { Session } from 'next-auth'
import { useRuntimeConfig } from '@/app/fragments/config'
import usePublicUrl from '@/utils/next-public-url'
import ThemeToggle from '@/app/components/ThemeToggle'

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

export default function Layout({ children }: React.PropsWithChildren) {
  const publicUrl = usePublicUrl()
  const { data: session } = Auth.useSession()
  const runtimeConfig = useRuntimeConfig()
  return (
    <div className="drawer min-w-screen min-h-screen">
      <input id="my-drawer-3" type="checkbox" className="drawer-toggle" />
      <div className="drawer-content flex flex-col">
        <div className="w-full navbar bg-primary gap-2">
          <div className="flex-grow flex-none">
            <div className="flex-none lg:hidden">
              <label htmlFor="my-drawer-3" className="btn btn-square btn-ghost">
                <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" className="inline-block w-6 h-6 stroke-black dark:stroke-white"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M4 6h16M4 12h16M4 18h16"></path></svg>
              </label>
            </div>
            <Link href={runtimeConfig.NEXT_PUBLIC_LANDING_PAGE}>
              <h1 className="text-black dark:text-white text-4xl font-bold p-2 cursor-pointer">P<span className="text-2xl">laybook</span> W<span className="text-2xl">orkflow</span> B<span className="text-2xl">uilder</span></h1>
            </Link>
          </div>
          <div className="navbar-end hidden md:flex">
            <div className="hidden lg:flex">
              <Link href="/playbooks"><button className="btn btn-ghost text-black dark:text-white">Published Playbooks</button></Link>
              <Link href="/explore"><button className="btn btn-ghost text-black dark:text-white">Explore Components</button></Link>
              <a className="text-black hover:text-black dark:text-white dark:hover:text-white" href="https://github.com/nih-cfde/playbook-partnership/blob/main/docs/user/index.md" target="_blank"><button className="btn btn-ghost">User Guide</button></a>
            </div>
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
              : <Link href="/api/auth/signin"><button className="btn btn-ghost text-black hover:text-black dark:text-white dark:hover:text-white">Sign in</button></Link>}
          </div>
        </div>

        {children}

        <div className="grid grid-flow-row md:grid-cols-3 sm:grid-cols-2 grid-cols-1 gap-4 justify-items-center items-center text-center mb-2">
          <div className="flex flex-col grid-cols-1">
            <a className="prose prose-sm" href="mailto:avi.maayan@mssm.edu">Contact Us</a>
            <a className="prose prose-sm" href="https://github.com/nih-cfde/playbook-partnership/blob/main/LICENSE" target="_blank">Usage License</a>
            <a className="prose prose-sm" href="https://github.com/nih-cfde/playbook-partnership/blob/main/docs/contributions.md" target="_blank">Contribute</a>
          </div>
          <div className="grid-cols-1">
            <a href="https://www.nih-cfde.org/" target="_blank">
              <img className="rounded h-20 dark:bg-white" src={`${publicUrl}/logos/CFDE.png`} />
            </a>
          </div>
          <div className="flex flex-col grid-cols-1 gap-1">
            <a className="btn btn-xs btn-secondary rounded-lg gap-1" href="https://github.com/nih-cfde/playbook-partnership" target="_blank">
              <img className="rounded-md w-4 justify-self-start" src={`${publicUrl}/GitHub-Mark.png`} />
              <span className="flex-grow text-black dark:text-white">View source code</span>
            </a>
            <a className="btn btn-xs btn-secondary rounded-lg gap-1" href="https://github.com/nih-cfde/playbook-partnership/issues/new" target="_blank">
              <img className="rounded-md w-4 justify-self-start" src={`${publicUrl}/GitHub-Mark.png`} />
              <span className="flex-grow text-black dark:text-white">Submit an issue</span>
            </a>
            <div className="flex self-center gap-2 items-center">
              <span className="prose">Light</span>
              <ThemeToggle />
              <span className="prose">Dark</span>
            </div>
          </div>
        </div>
      </div>
      <div className="drawer-side">
        <label htmlFor="my-drawer-3" className="drawer-overlay"></label>
        <ul className="menu p-4 w-80 bg-base-100">
          <li className="menu-title">
            <span>Playbook</span>
          </li>
          <li><Link href="/playbooks">Published Playbooks</Link></li>
          <li><Link href="/explore">Explore Components</Link></li>
          <li><a href="https://github.com/nih-cfde/playbook-partnership/blob/main/docs/user/index.md" target="_blank">User Guide</a></li>
          <li className="menu-title">
            <span>Account</span>
          </li>
          {session && session.user ?
              <>
                <li><Link href="/account">Settings</Link></li>
                <li><Link href="/api/auth/signout">Sign Out</Link></li>
              </>
              : <li><Link href="/api/auth/signin">Sign in</Link></li>}
        </ul>
      </div>
    </div>
  )
}