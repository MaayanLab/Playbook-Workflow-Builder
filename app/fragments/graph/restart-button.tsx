import React from 'react'
import { useRouter } from 'next/router'
import dynamic from 'next/dynamic'
import { start_icon, restart_icon } from '@/icons'

const Icon = dynamic(() => import('@/app/components/icon'))
const Bp4Alert = dynamic(() => import('@blueprintjs/core').then(({ Alert }) => Alert))

export default function RestartButton({ session_id }: { session_id?: string }) {
  const router = useRouter()
  const [isOpen, setIsOpen] = React.useState(false)
  const onConfirm = React.useCallback(() => {
    setIsOpen(false)
    router.push(`${session_id ? `/session/${session_id}` : ''}/graph/extend`, undefined, { shallow: true })
  }, [router])
  const disabled = router.asPath.endsWith('/graph/start') || router.asPath.endsWith('/graph/extend') || router.asPath.endsWith('/graph/start/extend')
  return (
    <>
      <button
        className="bp4-button bp4-minimal"
        disabled={disabled}
        onClick={evt => {
          if (evt.shiftKey) {
            onConfirm()
          } else {
            setIsOpen(true)
          }
        }}
      >
        <Icon icon={restart_icon} className={disabled ? 'fill-gray-400' : 'fill-black dark:fill-white'} />
      </button>
      <Bp4Alert
        cancelButtonText="Cancel"
        confirmButtonText="Restart"
        icon="reset"
        intent="warning"
        isOpen={isOpen}
        canEscapeKeyCancel
        canOutsideClickCancel
        onCancel={() => {setIsOpen(false)}}
        onConfirm={() => {onConfirm()}}
      >
        <p className="prose">
          Are you sure you want to restart? If you haven't saved the this session, it might be deleted.
        </p>
        <p className="prose">
          If you want to continue this session with additional inputs, you should instead click the <Icon icon={start_icon} /> icon in the breadcrumbs.
        </p>
        <p className="prose">
          <b>Tip:</b> Hold shift when clicking to skip this confirmation.
        </p>
      </Bp4Alert>
    </>
  )
}