import dynamic from 'next/dynamic'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import { save_icon } from '@/icons'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function SaveButton({ userPlaybook, toggleSave, updateRequired }: {
  userPlaybook?: { public: boolean },
  toggleSave: () => void,
  updateRequired: boolean,
}) {
  const session = useSessionWithId()
  const disabled = session.data?.user === undefined
  return (
    <>
      <button
        className="bp4-button bp4-minimal"
        disabled={disabled}
        onClick={() => { toggleSave() }}
      >
        <Icon
          icon={save_icon}
          color={disabled ? 'grey' : !userPlaybook ? 'black' : updateRequired ? 'crimson' : 'green'}
          title={disabled ? 'Sign In to Save' : 'Save'}
        />
      </button>
    </>
  )
}
