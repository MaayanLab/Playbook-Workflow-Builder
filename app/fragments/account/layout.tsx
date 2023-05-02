import * as Auth from 'next-auth/react'
import { useRouter } from 'next/router'
import React from 'react'
import { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'
import dynamic from 'next/dynamic'
import { Tabs as Bp4Tabs, Tab as Bp4Tab } from '@blueprintjs/core'

const Bp4Icon = dynamic(() => import('@blueprintjs/core').then(({ Icon }) => Icon))

const AccountProfile = dynamic(() => import('@/app/fragments/account/profile'))
const AccountSettings = dynamic(() => import('@/app/fragments/account/settings'))
const AccountUploads = dynamic(() => import('@/app/fragments/account/uploads'))
const AccountPlaybooks = dynamic(() => import('@/app/fragments/account/playbooks'))
const AccountSuggestions = dynamic(() => import('@/app/fragments/account/suggestions'))
const AccountBioCompute = dynamic(() => import('@/app/fragments/account/biocompute'))
const AccountCAVATICA = dynamic(() => import('@/app/fragments/account/cavatica'))

export default function Layout({ session }: { session: SessionWithId }) {
  const router = useRouter()
  const tab = (router.query.tab as string) || 'profile'
  return (
    <Bp4Tabs
      id="account-tabs"
      className="flex-grow"
      selectedTabId={tab}
      onChange={(tab) => {
        if (tab === 'signout') {
          Auth.signOut({ callbackUrl: '/account' })
        } else {
          router.push(`/account/${tab}`, undefined, { shallow: true })
        }
      }}
      renderActiveTabPanelOnly
      vertical
    >
      <span className='font-bold'>{session.user.email}</span>
      <Bp4Tab id="profile" title={<><Bp4Icon icon="person" /> Profile</>} panelClassName="flex-grow flex flex-col" panel={<AccountProfile session={session}  />} />
      <Bp4Tab id="settings" title={<><Bp4Icon icon="cog" /> Settings</>} panelClassName="flex-grow flex flex-col" panel={<AccountSettings session={session}  />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Data</span>
      <Bp4Tab id="uploads" title={<><Bp4Icon icon="upload" /> Uploads</>} panelClassName="flex-grow flex flex-col" panel={<AccountUploads />} />
      <Bp4Tab id="playbooks" title={<><Bp4Icon icon="control" /> Playbooks</>} panelClassName="flex-grow flex flex-col" panel={<AccountPlaybooks />} />
      <Bp4Tab id="suggestions" title={<><Bp4Icon icon="lightbulb" /> Suggestions</>} panelClassName="flex-grow flex flex-col" panel={<AccountSuggestions />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Integrations</span>
      <Bp4Tab id="biocompute" title={<><Bp4Icon icon="application" /> BioCompute</>} panelClassName="flex-grow flex flex-col" panel={<AccountBioCompute />} />
      <Bp4Tab id="cavatica" title={<><Bp4Icon icon="application" /> CAVATICA</>} panelClassName="flex-grow flex flex-col" panel={<AccountCAVATICA />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Session</span>
      <Bp4Tab id="signout" title={<><Bp4Icon icon="log-out" /> Sign Out</>} />
    </Bp4Tabs>
  )
}
