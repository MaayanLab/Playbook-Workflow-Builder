import React from 'react'
import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query'
import { SessionWithId } from '@/app/pages/api/auth/[...nextauth]'
import { useTRPC } from '@/app/trpc/client'
import dynamic from 'next/dynamic'
import classNames from 'classnames'

const Bp5ControlGroup = dynamic(() => import('@blueprintjs/core').then(({ ControlGroup }) => ControlGroup))
const Bp5FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const Bp5InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))
const Bp5Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export default function Profile({ session }: { session: SessionWithId }) {
  const trpc = useTRPC()
  const queryClient = useQueryClient()
  const { data: userProfile, isLoading } = useQuery(trpc.users.me.queryOptions())
  const updateProfile = useMutation(trpc.users.update.mutationOptions({
    onSuccess: () => queryClient.invalidateQueries({ queryKey: trpc.users.me.queryKey() }),
  }))
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
      <h3 className="bp5-heading">Profile Settings</h3>
      <progress className={classNames('progress w-full', { 'hidden': !(isLoading || updateProfile.isPending) })}></progress>
      <form onSubmit={async (evt) => {
        evt.preventDefault()
        updateProfile.mutate(userProfileDraft)
      }} method="POST">
        <Bp5FormGroup
          label="Authorship Information"
          helperText="Let us know who you are and how to contact you"
        >
          <Bp5ControlGroup fill vertical>
            <Bp5ControlGroup fill>
              <Bp5InputGroup
                type="text"
                placeholder="Name"
                value={userProfileDraft.name || ''}
                onChange={evt => {
                  setUserProfileDraft(({ ...user }) => ({ ...user, name: evt.target.value }))
                }}
                leftIcon="person"
              />
              <Bp5InputGroup
                type="email"
                placeholder="Email"
                readOnly
                value={session.user?.email || ''}
                leftIcon="envelope"
              />
            </Bp5ControlGroup>
            <Bp5InputGroup
              type="text"
              placeholder="Affiliation"
              value={userProfileDraft.affiliation || ''}
              onChange={evt => {
                setUserProfileDraft(({ ...user }) => ({ ...user, affiliation: evt.target.value }))
              }}
              leftIcon="office"
            />
          </Bp5ControlGroup>
        </Bp5FormGroup>
        <Bp5Button
          type="submit"
          intent="success"
        >Update Profile</Bp5Button>
      </form>
    </>
  )
}
