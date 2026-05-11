import React from 'react'
import { ExLink } from '@/app/fragments/ex-router'
import * as Auth from 'next-auth/react'
import { useRuntimeConfig } from '@/app/fragments/config'
import usePublicUrl from '@/utils/next-public-url'
import ThemeToggle from '@/app/components/ThemeToggle'
import { Waypoint, Waypoints } from '@/app/components/waypoint'
import dynamic from 'next/dynamic'
import classNames from 'classnames'

const UserAvatar = dynamic(() => import('@/app/fragments/playbook/avatar'))
const DismissableAlert = dynamic(() => import('@/app/fragments/DismissableAlert'))

export default function Layout({ children, sidebar }: React.PropsWithChildren<{ sidebar?: React.ReactNode }>) {
  const publicUrl = usePublicUrl()
  const { data: session } = Auth.useSession()
  const runtimeConfig = useRuntimeConfig()
  const [isSidebarCollapsed, setIsSidebarCollapsed] = React.useState(false)
  return (
    <div className="drawer w-screen h-screen">
      <input id="my-drawer-3" type="checkbox" className="drawer-toggle" />
      <Waypoints className="drawer-content flex flex-col w-screen h-screen overflow-hidden">
        <Waypoint id="top" />
        <div className="navbar bg-primary gap-2">
          <div className="flex-grow flex-none">
            <div className="flex-none xl:hidden">
              <label htmlFor="my-drawer-3" className="btn btn-square btn-ghost">
                <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" className="inline-block w-6 h-6 stroke-black dark:stroke-white"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M4 6h16M4 12h16M4 18h16"></path></svg>
              </label>
            </div>
            <ExLink
              href={runtimeConfig.NEXT_PUBLIC_LANDING_PAGE}
              className="flex flex-row items-center"
            >
              <img
                className="p-2 max-w-16 cursor-pointer dark:stroke-white"
                src={`${publicUrl}/PWB-logo.svg`}
              />
              <h1 className="text-black dark:text-white text-lg sm:text-2xl md:text-4xl font-bold p-2 cursor-pointer"><span className="text-[140%]">P</span>laybook <span className="text-[140%]">W</span>orkflow <span className="text-[140%]">B</span>uilder</h1>
            </ExLink>
          </div>
          <div className="navbar-end hidden md:flex">
            <div className="hidden xl:flex">
              <ExLink href="/playbooks"><button className="btn btn-ghost text-black dark:text-white">Published Playbooks</button></ExLink>
              <ExLink href="/chat"><button className="btn btn-ghost text-black dark:text-white">Text to Workflow</button></ExLink>
              <ExLink href="/components"><button className="btn btn-ghost text-black dark:text-white">Component Catalog</button></ExLink>
              <ExLink href="/explore"><button className="btn btn-ghost text-black dark:text-white">Component Graph</button></ExLink>
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
                  <li className="menu-title prose"><span>Data</span></li>
                  <li><ExLink href="/account/uploads">Uploads</ExLink></li>
                  <li><ExLink href="/account/playbooks">Playbooks</ExLink></li>
                  <li><ExLink href="/account/suggestions">Suggestions</ExLink></li>
                  <li className="menu-title prose"><span>Account</span></li>
                  <li><ExLink href="/account">Settings</ExLink></li>
                  <li><ExLink href="/api/auth/signout" onClick={() => {Auth.signOut()}}>Sign Out</ExLink></li>
                </ul>
              </div>
              : <button className="btn btn-ghost text-black hover:text-black dark:text-white dark:hover:text-white" onClick={() => {Auth.signIn()}}>Sign in</button>}
          </div>
        </div>

        <div className="grow flex flex-col sm:flex-row-reverse justify-stretch overflow-hidden">
          <div className="grow flex flex-col sm:h-full overflow-auto">
            <div className="container mx-auto my-2 flex flex-col gap-2">
              <DismissableAlert id="plos" style={{ backgroundColor: '#d4fb79' }}>
                <ExLink href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012901" className="dark:text-inherit dark:hover:text-inherit">
                  An article that describes the Playbook Workflow Builder project is now published in PLOS Computational Biology!
                </ExLink>
              </DismissableAlert>
            </div>
            
            {children}

            <div className="grid grid-flow-row md:grid-cols-3 sm:grid-cols-2 grid-cols-1 gap-4 justify-items-center items-center text-center mb-2">
              <div className="flex flex-col grid-cols-1">
                <a className="prose prose-sm max-w-none" href="mailto:avi.maayan@mssm.edu">Contact Us</a>
                <a className="prose prose-sm max-w-none" href="https://github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/LICENSE" target="_blank">Usage License</a>
                <a className="prose prose-sm max-w-none" href="https://github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/docs/contributions.md" target="_blank">Contribute</a>
              </div>
              <div className="grid-cols-1">
                <a href="https://info.cfde.cloud" target="_blank">
                  <img className="rounded h-20 dark:bg-white" src={`${publicUrl}/logos/CFDE.png`} />
                </a>
              </div>
              <div className="flex flex-col grid-cols-1 gap-1">
                <a className="btn btn-xs btn-secondary rounded-lg gap-1" href="https://github.com/MaayanLab/Playbook-Workflow-Builder" target="_blank">
                  <img className="rounded-md w-4 justify-self-start" src={`${publicUrl}/GitHub-Mark.png`} />
                  <span className="flex-grow text-black dark:text-white">View source code</span>
                </a>
                <a className="btn btn-xs btn-secondary rounded-lg gap-1" href="https://github.com/MaayanLab/Playbook-Workflow-Builder/issues/new" target="_blank">
                  <img className="rounded-md w-4 justify-self-start" src={`${publicUrl}/GitHub-Mark.png`} />
                  <span className="flex-grow text-black dark:text-white">Submit an issue</span>
                </a>
                <div className="flex self-center gap-2 items-center">
                  <span className="prose max-w-none">Light</span>
                  <ThemeToggle />
                  <span className="prose max-w-none">Dark</span>
                </div>
              </div>
            </div>
            <div className="prose prose-sm max-w-none prose-p:m-0 alert bg-secondary flex-col gap-1 rounded-none">
              <p className="flex-col font-bold">Please acknowledge Playbook Workflow Builder in your publications by citing the following reference:</p>
              <p className="inline text-center">Clarke, D. J. B. et al. Playbook workflow builder: Interactive construction of bioinformatics workflows. PLOS Computational Biology vol. 21 e1012901 (2025). <a href="https://doi.org/10.1371/journal.pcbi.1012901" target="_blank">doi:10.1371/journal.pcbi.1012901</a></p>
            </div>
            <Waypoint id="bottom" />
          </div>
          {sidebar && <>
            <div className="flex flex-col sm:flex-row-reverse">
              <div className={classNames("self-center relative")}>
                <button className={classNames("absolute btn btn-xs sm:rotate-90 sm:active:focus:rotate-90 whitespace-nowrap sm:-left-8 left-[25vw]", { "-top-4 sm:top-0 sm:-left-9": !isSidebarCollapsed,  "-top-6 sm:top-0 sm:-left-2": isSidebarCollapsed })}  onClick={evt => {setIsSidebarCollapsed(!isSidebarCollapsed)}}>
                  {isSidebarCollapsed ? 'Chat' : 'Hide Chat'}
                </button>
              </div>
              <div className={classNames("w-full h-[50vh] sm:h-full sm:w-[max(33vw,20em)] shrink-0 border-t sm:border-t-0 sm:border-r flex flex-col overflow-hidden", { 'hidden': isSidebarCollapsed })}>{sidebar}</div>
            </div>
          </>}
        </div>

      </Waypoints>
      <div className="drawer-side">
        <label htmlFor="my-drawer-3" className="drawer-overlay"></label>
        <ul className="menu p-4 w-80 bg-base-100">
          <li className="menu-title prose">
            <span>Playbook</span>
          </li>
          <li><ExLink href="/playbooks">Published Playbooks</ExLink></li>
          <li><ExLink href="/chat">Text to Workflow</ExLink></li>
          <li><ExLink href="/components">Component Catalog</ExLink></li>
          <li><ExLink href="/explore">Component Graph</ExLink></li>
          <li><a href="https://github.com/MaayanLab/Playbook-Workflow-Builder/blob/main/docs/user/index.md" target="_blank">User Guide</a></li>
          {session && session.user ?
              <>
                    <li className="menu-title prose">
                      <span>Data</span>
                    </li>
                    <li><ExLink href="/account/uploads">Uploads</ExLink></li>
                    <li><ExLink href="/account/playbooks">Playbooks</ExLink></li>
                    <li><ExLink href="/account/suggestions">Suggestions</ExLink></li>
                    <li className="menu-title prose">
                      <span>Account</span>
                    </li>
                    <li><ExLink href="/account">Settings</ExLink></li>
                    <li><ExLink href="/api/auth/signout" onClick={() => {Auth.signOut()}}>Sign Out</ExLink></li>
              </>
              :
              <>
                <li className="menu-title prose">
                  <span>Account</span>
                </li>
                <li><ExLink href="/api/auth/signin" onClick={() => {Auth.signIn()}}>Sign in</ExLink></li>
              </>}
        </ul>
      </div>
    </div>
  )
}