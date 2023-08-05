import React from 'react'
import useSWR from 'swr'
import useSWRMutation from 'swr/mutation'
import { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'
import fetcher from '@/utils/next-rest-fetcher'
import dynamic from 'next/dynamic'
import classNames from 'classnames'

const Bp4ControlGroup = dynamic(() => import('@blueprintjs/core').then(({ ControlGroup }) => ControlGroup))
const Bp4FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const Bp4InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))
const Bp4Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

const poster = (endpoint: string, { arg }: { arg: any }) => fetch(endpoint, { method: 'POST', body: JSON.stringify(arg) }).then(res => res.json())

export default function Profile({ session }: { session: SessionWithId }) {
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
      <progress className={classNames('progress w-full', { 'hidden': !(isLoading || isMutating) })}></progress>
      <form onSubmit={async (evt) => {
        evt.preventDefault()
        setUserProfile(userProfileDraft, { revalidate: false })
      }} method="POST">
        <Bp4FormGroup
          label="Authorship Information"
          helperText="Let us know who you are and how to contact you"
        >
          <Bp4ControlGroup fill vertical>
            <Bp4ControlGroup fill>
              <Bp4InputGroup
                type="text"
                placeholder="Name"
                value={userProfileDraft.name || ''}
                onChange={evt => {
                  setUserProfileDraft(({ ...user }) => ({ ...user, name: evt.target.value }))
                }}
                leftIcon="person"
              />
              <Bp4InputGroup
                type="email"
                placeholder="Email"
                readOnly
                value={session.user?.email || ''}
                leftIcon="envelope"
              />
            </Bp4ControlGroup>
            <Bp4InputGroup
              type="text"
              placeholder="Affiliation"
              value={userProfileDraft.affiliation || ''}
              onChange={evt => {
                setUserProfileDraft(({ ...user }) => ({ ...user, affiliation: evt.target.value }))
              }}
              leftIcon="office"
            />
          </Bp4ControlGroup>
        </Bp4FormGroup>
        <Bp4Button
          type="submit"
          intent="success"
        >Update Profile</Bp4Button>
      </form>
    </>
  )
}
