import dynamic from 'next/dynamic'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import { save_icon } from '@/icons'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function SaveButton({ userPlaybook, toggleSave, updateRequired }: {
  userPlaybook?: { public: boolean },
  toggleSave: () => void,
  updateRequired: boolean,
}) {
  const userSession = useSessionWithId()
  const disabled = !userSession.data?.user
  return (
    <>
      <button
        className="bp4-button bp4-minimal"
        disabled={disabled}
        onClick={() => { toggleSave() }}
      >
        <Icon
          icon={save_icon}
          className={disabled ? 'fill-gray-400' : !userPlaybook ? 'fill-black dark:fill-white' : updateRequired ? 'fill-red-600' : 'fill-green-600'}
          title={disabled ? 'Sign In to Save' : 'Save'}
        />
      </button>
    </>
  )
}
