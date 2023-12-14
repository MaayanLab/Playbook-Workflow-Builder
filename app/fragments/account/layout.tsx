import * as Auth from 'next-auth/react'
import { useRouter } from 'next/router'
import React from 'react'
import { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'
import dynamic from 'next/dynamic'
import { Tabs as Bp5Tabs, Tab as Bp5Tab } from '@blueprintjs/core'

const Bp5Icon = dynamic(() => import('@blueprintjs/core').then(({ Icon }) => Icon))

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
    <Bp5Tabs
      id="account-tabs"
      className="flex-grow flex overflow-hidden"
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
      <span className='font-bold'>{session.user?.email ?? null}</span>
      <Bp5Tab id="profile" title={<><Bp5Icon icon="person" /> Profile</>} panelClassName="flex-grow flex flex-col overflow-hidden" panel={<AccountProfile session={session}  />} />
      <Bp5Tab id="settings" title={<><Bp5Icon icon="cog" /> Settings</>} panelClassName="flex-grow flex flex-col overflow-hidden" panel={<AccountSettings session={session}  />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Data</span>
      <Bp5Tab id="uploads" title={<><Bp5Icon icon="upload" /> Uploads</>} panelClassName="flex-grow flex flex-col overflow-hidden" panel={<AccountUploads />} />
      <Bp5Tab id="playbooks" title={<><Bp5Icon icon="control" /> Playbooks</>} panelClassName="flex-grow flex flex-col overflow-hidden" panel={<AccountPlaybooks />} />
      <Bp5Tab id="suggestions" title={<><Bp5Icon icon="lightbulb" /> Suggestions</>} panelClassName="flex-grow flex flex-col overflow-hidden" panel={<AccountSuggestions />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Integrations</span>
      <Bp5Tab id="biocompute" title={<><Bp5Icon icon="application" /> BioCompute</>} panelClassName="flex-grow flex flex-col overflow-hidden" panel={<AccountBioCompute />} />
      <Bp5Tab id="cavatica" title={<><Bp5Icon icon="application" /> CAVATICA</>} panelClassName="flex-grow flex flex-col overflow-hidden" panel={<AccountCAVATICA />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Session</span>
      <Bp5Tab id="signout" title={<><Bp5Icon icon="log-out" /> Sign Out</>} />
    </Bp5Tabs>
  )
}
