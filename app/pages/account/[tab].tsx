import React from 'react'
import Head from 'next/head'
import * as Auth from 'next-auth/react'
import dynamic from 'next/dynamic'
import { Session } from 'next-auth'
import { Alert, ControlGroup, FormGroup, Icon, InputGroup, Tab, Tabs } from '@blueprintjs/core'

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
        {(session && session.user) ?
          <AccountUI session={session as SessionWithUser} />
          : <AccountSignIn />}
      </main>

      <Footer />
    </div>
  )
}

function AccountSignIn() {
  return (
    <div className="flex-grow flex items-center justify-center">
      <div className="text-center">
        <p className="mb-2">You must sign in to view this page.</p>
        <Button onClick={() => Auth.signIn()}>Sign in</Button>
      </div>
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
      <Tab id="profile" title={<><Icon icon="person" /> Profile</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIProfile session={session}  />} />
      <Tab id="account" title={<><Icon icon="cog" /> Account</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIAccount session={session}  />} />
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

function AccountUIProfile({ session }: { session: SessionWithUser }) {
  const [userDraft, setUserDraft] = React.useState(session.user)
  return (
    <>
      <h3 className="bp4-heading">Profile Settings</h3>
      <div>
        <FormGroup
          label="Authorship Information"
          helperText="Let us know who you are and how to contact you"
        >
          <ControlGroup fill vertical>
            <ControlGroup fill>
              <InputGroup
                type="text"
                placeholder="Name"
                value={userDraft.name || ''}
                onChange={evt => {
                  setUserDraft(({ ...user }) => ({ ...user, name: evt.target.value }))
                }}
                leftIcon="person"
              />
              <InputGroup
                type="email"
                placeholder="Email"
                readOnly
                value={userDraft.email || ''}
                onChange={evt => {
                  setUserDraft(({ ...user }) => ({ ...user, email: evt.target.value }))
                }}
                leftIcon="envelope"
              />
            </ControlGroup>
            <InputGroup
              type="text"
              placeholder="Affiliation"
              value={userDraft.org || ''}
              onChange={evt => {
                setUserDraft(({ ...user }) => ({ ...user, org: evt.target.value }))
              }}
              leftIcon="office"
            />
          </ControlGroup>
        </FormGroup>
        <Button
          intent="success"
          onClick={() => {
            // TODO: Update profile
          }}
        >Update Profile</Button>
      </div>
    </>
  )
}

function AccountUIAccount({ session }: { session: SessionWithUser }) {
  return (
    <>
      <div className="mb-2"><AccountUIAccountDeleteAccount session={session} /></div>
    </>
  )
}

function AccountUIAccountDeleteAccount({ session }: { session: SessionWithUser }) {
  const [deletionConfirmation, setDeletionConfirmation] = React.useState(false)
  return (
    <>
      <h3 className="bp4-heading text-red-600">Delete Account</h3>
      <Button
        intent="danger"
        onClick={() => setDeletionConfirmation(true)}
      >Delete your account</Button>
      <Alert
        cancelButtonText="Cancel"
        confirmButtonText="Delete Account"
        icon="delete"
        intent="danger"
        isOpen={deletionConfirmation}
        canEscapeKeyCancel
        canOutsideClickCancel
        onCancel={() => {setDeletionConfirmation(false)}}
        onConfirm={() => {
          // TODO: actually delete account
        }}
      >
        Are you sure you want to delete your account?
        After clicking Delete Account, your account <b>and all associated content including uploads and saved graphs</b> will be subject to deletion and <b>cannot be restored</b>.<br />
      </Alert>
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
