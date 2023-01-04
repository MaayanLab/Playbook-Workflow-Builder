import React from 'react'
import Head from 'next/head'
import * as Auth from 'next-auth/react'
import dynamic from 'next/dynamic'
import { Session } from 'next-auth'
import { Alert, ControlGroup, FormGroup, Icon, InputGroup, Tab, Tabs } from '@blueprintjs/core'
import { useRouter } from 'next/router'
import Image from 'next/image'

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
  const router = useRouter()
  const tab = (router.query.tab as string) || 'profile'
  return (
    <Tabs
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
      <Tab id="profile" title={<><Icon icon="person" /> Profile</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIProfile session={session}  />} />
      <Tab id="settings" title={<><Icon icon="cog" /> Settings</>} panelClassName="flex-grow flex flex-col" panel={<AccountUISettings session={session}  />} />
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

function AccountUISettings({ session }: { session: SessionWithUser }) {
  return (
    <>
      <div className="mb-2"><AccountUISettingsDeleteAccount session={session} /></div>
    </>
  )
}

function AccountUISettingsDeleteAccount({ session }: { session: SessionWithUser }) {
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

import CAVATICAProjectGuide from '@/app/public/CAVATICA-guide-project.png'
import CAVATICAAPIKeyGuide from '@/app/public/CAVATICA-guide-apikey.png'

function AccountUICAVATICA() {
  return (
    <>
      <h3 className="bp4-heading">CAVATICA Integration</h3>
      <p className="bp4-running-text">
        <a href="https://www.cavatica.org/">CAVATICA</a> is a data analysis and sharing platform designed to accelerate discovery in a scalable, cloud-based compute environment where data, results, and workflows are shared among the world's research community. Developed by Seven Bridges and funded in-part by a grant from the National Institutes of Health (NIH) Common Fund, CAVATICA is continuously updated with new tools and datasets.</p>
        <p className="my-2">CAVATICA offers secure storage and compute in a cloud environment. Appyter integration with CAVATICA enables you to execute Appyters against data in a CAVATICA project and use CAVATICA-managed computational resources.</p>
        <p className="my-2">To use CAVATICA, you must <a href="https://pgc-accounts.sbgenomics.com/auth/login" target="_blank">login or create an account</a>, then a project should be established which will be used for file storage and executions that are configured.</p>
        <Image src={CAVATICAProjectGuide} alt="CAVATICA Project Creation Guide" />
        <p className="my-2">To use CAVATICA you must register your CAVATICA API Key with Appyters, this key can be located at <a href="https://cavatica.sbgenomics.com/developer/token" target="_blank">https://cavatica.sbgenomics.com/developer/token</a></p>
        <Image src={CAVATICAAPIKeyGuide} alt="CAVATICA API Key Guide" />
        <div className="my-4">
          <FormGroup
            label="API Key"
            labelFor="cavatica-api-key"
            helperText="Provide your CAVATICA credentials"
          >
            <InputGroup
              id="cavatica-api-key"
              type="string"
              placeholder="e.g. 08cd35123..."
              leftIcon="key"
            />
          </FormGroup>
          <FormGroup
            label="CAVATICA Project"
            labelFor="cavatica-project"
            helperText="Specify the default project in CAVATICA"
          >
            <InputGroup
              id="cavatica-project"
              type="string"
              placeholder="e.g. youruser/yourproject"
              leftIcon="projects"
            />
          </FormGroup>
          <ControlGroup>
            <Button
              intent="success"
              onClick={() => {
                // TODO
              }}>Save</Button>
            <Button
              intent="danger"
              onClick={() => {
                // TODO
              }}>Delete</Button>
          </ControlGroup>
        </div>
    </>
  )
}
