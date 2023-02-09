import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import CAVATICAAPIKeyGuide from '@/app/public/CAVATICA-guide-apikey.png'
import CAVATICAProjectGuide from '@/app/public/CAVATICA-guide-project.png'
import { Alert, ControlGroup, FormGroup, Icon as Bp4Icon, InputGroup, Tab, Tabs } from '@blueprintjs/core'
import * as Auth from 'next-auth/react'
import * as schema from '@/db'
import type { TypedSchemaRecord } from '@/spec/sql'
import dynamic from 'next/dynamic'
import Head from 'next/head'
import Image from 'next/image'
import { useRouter } from 'next/router'
import React from 'react'
import useSWR from 'swr'
import useSWRMutation from 'swr/mutation'
import { SessionWithId } from '../api/auth/[...nextauth]'
import Icon from '@/app/components/icon'
import { delete_icon, fork_icon } from '@/icons'
import { z } from 'zod'
import { FileInput, FileURL } from '@/components/core/file'

const Header = dynamic(() => import('@/app/fragments/playbook/header'))
const Footer = dynamic(() => import('@/app/fragments/playbook/footer'))
const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

const fetcher = (endpoint: string) => fetch(endpoint).then(res => res.json())
const poster = (endpoint: string, { arg }: { arg: any }) => fetch(endpoint, { method: 'POST', body: JSON.stringify(arg) }).then(res => res.json())
const deleter = (endpoint: string, { arg }: { arg: any }) => fetch(`${endpoint}/${arg}/delete`, { method: 'POST' }).then(res => res.json())

export default function Account() {
  const { data: session } = useSessionWithId({ required: true })
  return (
    <div className="flex flex-col min-w-screen min-h-screen">
      <Head>
        <title>Playbook: Account</title>
      </Head>

      <Header />
      
      <main className="flex-grow container mx-auto py-4 flex flex-row">
        {(session && session.user) ?
          <AccountUI session={session} />
          : null}
      </main>

      <Footer />
    </div>
  )
}

function AccountUI({ session }: { session: SessionWithId }) {
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
      <Tab id="profile" title={<><Bp4Icon icon="person" /> Profile</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIProfile session={session}  />} />
      <Tab id="settings" title={<><Bp4Icon icon="cog" /> Settings</>} panelClassName="flex-grow flex flex-col" panel={<AccountUISettings session={session}  />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Data</span>
      <Tab id="uploads" title={<><Bp4Icon icon="upload" /> Uploads</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIUploads />} />
      <Tab id="playbooks" title={<><Bp4Icon icon="control" /> Playbooks</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIPlaybooks />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Integrations</span>
      <Tab id="biocompute" title={<><Bp4Icon icon="application" /> BioCompute</>} panelClassName="flex-grow flex flex-col" panel={<AccountUIBioCompute />} />
      <Tab id="cavatica" title={<><Bp4Icon icon="application" /> CAVATICA</>} panelClassName="flex-grow flex flex-col" panel={<AccountUICAVATICA />} />
      <hr className="h-px my-1 border-0 bg-secondary w-full" />
      <span className='font-bold'>Session</span>
      <Tab id="signout" title={<><Bp4Icon icon="log-out" /> Sign Out</>} />
    </Tabs>
  )
}

function AccountUIProfile({ session }: { session: SessionWithId }) {
  const { data: userProfile, isLoading } = useSWR<{ name: string, affiliation: string }>(`/api/db/user/profile`, fetcher)
  const { trigger: setUserProfile, isMutating } = useSWRMutation('/api/db/user/profile', poster)
  const [userProfileDraft, setUserProfileDraft] = React.useState({
    name: '',
    affiliation: '',
  })
  React.useEffect(() => {
    if (userProfile) {
      setUserProfileDraft({ name: userProfile.name || '', affiliation: userProfile.affiliation || '' })
    }
  }, [userProfile])
  return (
    <>
      <h3 className="bp4-heading">Profile Settings</h3>
      <progress className={`progress w-full ${isLoading || isMutating ? '' : 'hidden'}`}></progress>
      <form onSubmit={async (evt) => {
        evt.preventDefault()
        setUserProfile(userProfileDraft, { revalidate: false })
      }} method="POST">
        <FormGroup
          label="Authorship Information"
          helperText="Let us know who you are and how to contact you"
        >
          <ControlGroup fill vertical>
            <ControlGroup fill>
              <InputGroup
                type="text"
                placeholder="Name"
                value={userProfileDraft.name || ''}
                onChange={evt => {
                  setUserProfileDraft(({ ...user }) => ({ ...user, name: evt.target.value }))
                }}
                leftIcon="person"
              />
              <InputGroup
                type="email"
                placeholder="Email"
                readOnly
                value={session.user.email || ''}
                leftIcon="envelope"
              />
            </ControlGroup>
            <InputGroup
              type="text"
              placeholder="Affiliation"
              value={userProfileDraft.affiliation || ''}
              onChange={evt => {
                setUserProfileDraft(({ ...user }) => ({ ...user, affiliation: evt.target.value }))
              }}
              leftIcon="office"
            />
          </ControlGroup>
        </FormGroup>
        <Button
          type="submit"
          intent="success"
        >Update Profile</Button>
      </form>
    </>
  )
}

function AccountUISettings({ session }: { session: SessionWithId }) {
  return (
    <>
      <div className="mb-2"><AccountUISettingsDeleteAccount session={session} /></div>
    </>
  )
}

function AccountUISettingsDeleteAccount({ session }: { session: SessionWithId }) {
  const [deletionConfirmation, setDeletionConfirmation] = React.useState(false)
  const { trigger: deleteUser, isMutating } = useSWRMutation('/api/db/user/delete', poster)
  return (
    <>
      <h3 className="bp4-heading text-red-600">Delete Account</h3>
      <progress className={`progress w-full ${isMutating ? '' : 'hidden'}`}></progress>
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
          deleteUser(undefined, { revalidate: false })
            .then(() => Auth.signOut({ callbackUrl: '/' }))
        }}
      >
        Are you sure you want to delete your account?
        After clicking Delete Account, your account <b>and all associated content including uploads and saved graphs</b> will be subject to deletion and <b>cannot be restored</b>.<br />
      </Alert>
    </>
  )
}

function humanSize(size: number) {
  const units = ['B', 'KB', 'MB', 'GB', 'TB']
  let unit
  for (let i = 0; i < units.length-1; i++) {
    unit = units[i]
    if (size < 1000) {
      break
    } else {
      size /= 1000
    }
  }
  return `${size.toFixed(2)} ${unit}`
}

function AccountUIUploads() {
  const router = useRouter()
  const { data: uploads, isLoading } = useSWR<Array<TypedSchemaRecord<typeof schema.user_upload_complete>>>('/api/db/user/uploads', fetcher)
  const [uploadToDelete, setUploadToDelete] = React.useState<TypedSchemaRecord<typeof schema.user_upload_complete> | undefined>(undefined)
  const { trigger: deleteUpload, isMutating } = useSWRMutation('/api/db/user/uploads', deleter)
  return (
    <>
      <h3 className="bp4-heading">Uploads</h3>
      <progress className={`progress w-full ${isLoading || isMutating ? '' : 'hidden'}`}></progress>
      {uploads ? (
        <div className="overflow-x-auto">
          <table className="table table-compact w-full">
            <thead>
              <tr>
                <th>Filename</th>
                <th>URL</th> 
                <th>sha256</th>
                <th>Size</th>
                <th>Timestamp</th>
                <th>Actions</th>
              </tr>
            </thead> 
            <tbody>
              {uploads.map(upload => (
                <tr key={upload.id}>
                  <td>{upload.filename}</td>
                  <td>{upload.url}</td>
                  <td>{upload.sha256.slice(0, 5)}...{upload.sha256.slice(-5)}</td>
                  <td>{humanSize(upload.size)}</td>
                  <td>{upload.created}</td>
                  <td className="flex flex-row">
                    <button onClick={async () => {
                      const req = await fetch(`/api/db/fpl/start/extend`, {
                        method: 'POST',
                        body: JSON.stringify({
                          type: FileInput.spec,
                          inputs: {},
                          data: {
                            type: FileURL.spec,
                            value: FileURL.codec.encode(upload.url),
                          },
                        })
                      })
                      const res = z.string().parse(await req.json())
                      router.push(`/graph/${res}/extend`)
                    }}>
                      <Icon icon={fork_icon} color="black" />
                    </button>
                    <button onClick={() => {
                      setUploadToDelete(upload)
                    }}>
                      <Icon icon={delete_icon} color="black" />
                    </button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      ) : null}
      <Alert
        cancelButtonText="Cancel"
        confirmButtonText="Delete Upload"
        icon="delete"
        intent="danger"
        isOpen={uploadToDelete !== undefined}
        canEscapeKeyCancel
        canOutsideClickCancel
        onCancel={() => {setUploadToDelete(undefined)}}
        onConfirm={() => {
          if (!uploadToDelete) return
          deleteUpload(uploadToDelete.id, { revalidate: true })
            .then(() => setUploadToDelete(undefined))
        }}
      >
        Are you sure you want to delete {uploadToDelete?.filename} uploaded at {uploadToDelete?.created.toString()}?
        After clicking Delete Upload, your upload will be subject to deletion and <b>cannot be restored</b>.<br />
      </Alert>
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
