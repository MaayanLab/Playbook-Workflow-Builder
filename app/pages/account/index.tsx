import React from 'react'
import Head from 'next/head'
import * as Auth from 'next-auth/react'
import dynamic from 'next/dynamic'
import { Session } from 'next-auth'
import { Icon, Tab, Tabs } from '@blueprintjs/core'

const Header = dynamic(() => import('@/app/fragments/playbook/header'))
const Footer = dynamic(() => import('@/app/fragments/playbook/footer'))
const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

type SessionWithUser = Session & { user: NonNullable<Session['user']> }

export default function Account() {
  const { data: session } = Auth.useSession()
  return (
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook: Account</title>
      </Head>

      <Header homepage="/account" />
      
      <main className="flex-grow container mx-auto py-4 flex flex-row">
        {(session && session.user) ? <AccountUI session={session as SessionWithUser} /> : <>
          Not signed in <br/>
          <Button onClick={() => Auth.signIn()}>Sign in</Button>
        </>}
      </main>

      <Footer />
    </div>
  )
}

function AccountUI({ session }: { session: SessionWithUser }) {
  const [tab, setTab] = React.useState('profile')
  return (
    <Tabs
      className="flex-grow"
      selectedTabId={tab}
      onChange={(tab) => {
        if (tab === 'signout') {
          Auth.signOut()
        } else {
          setTab(tab.toString())
        }
      }}
      renderActiveTabPanelOnly
      vertical
    >
      <span className='font-bold'>{session.user.email}</span>
      <Tab id="profile" title={<><Icon icon="person" /> Profile</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIProfile />} />
      <Tab id="account" title={<><Icon icon="cog" /> Account</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIAccount />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Data</span>
      <Tab id="uploads" title={<><Icon icon="upload" /> Uploads</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIUploads />} />
      <Tab id="playbooks" title={<><Icon icon="control" /> Playbooks</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIPlaybooks />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Integrations</span>
      <Tab id="biocompute" title={<><Icon icon="application" /> BioCompute</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIBioCompute />} />
      <Tab id="cavatica" title={<><Icon icon="application" /> CAVATICA</>} panelClassName="flex-grow flex flex-col" panel={<AccountUICAVATICA />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Session</span>
      <Tab id="signout" title={<><Icon icon="log-out" /> Sign Out</>} />
    </Tabs>
  )
}

function AccountUIProfile() {
  return (
    <>
      <h3 className="bp4-heading">Profile Settings</h3>
      TODO
    </>
  )
}

function AccountUIAccount() {
  return (
    <>
      <h3 className="bp4-heading">Account Settings</h3>
      TODO
    </>
  )
}

function AccountUIUploads() {
  return (
    <>
      <h3 className="bp4-heading">Uploads</h3>
      TODO
    </>
  )
}

function AccountUIPlaybooks() {
  return (
    <>
      <h3 className="bp4-heading">Playbooks</h3>
      TODO
    </>
  )
}


function AccountUIBioCompute() {
  return (
    <>
      <h3 className="bp4-heading">BioCompute Integration</h3>
      TODO
    </>
  )
}

function AccountUICAVATICA() {
  return (
    <>
      <h3 className="bp4-heading">CAVATICA Integration</h3>
      TODO
    </>
  )
}
