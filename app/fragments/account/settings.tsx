import React from 'react'
import { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'
import * as Auth from 'next-auth/react'
import useSWRMutation from 'swr/mutation'
import dynamic from 'next/dynamic'
import classNames from 'classnames'

const Bp4Alert = dynamic(() => import('@blueprintjs/core').then(({ Alert }) => Alert))
const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

const poster = (endpoint: string, { arg }: { arg: any }) => fetch(endpoint, { method: 'POST', body: JSON.stringify(arg) }).then(res => res.json())

function DeleteAccount({ session }: { session: SessionWithId }) {
  const [deletionConfirmation, setDeletionConfirmation] = React.useState(false)
  const { trigger: deleteUser, isMutating } = useSWRMutation('/api/db/user/delete', poster)
  return (
    <>
      <h3 className="bp4-heading text-red-600">Delete Account</h3>
      <progress className={classNames('progress w-full', { 'hidden': !isMutating })}></progress>
      <Bp4Button
        intent="danger"
        onClick={() => setDeletionConfirmation(true)}
      >Delete your account</Bp4Button>
      <Bp4Alert
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
      </Bp4Alert>
    </>
  )
}

export default function Settings({ session }: { session: SessionWithId }) {
  return (
    <>
      <div className="mb-2"><DeleteAccount session={session} /></div>
    </>
  )
}
