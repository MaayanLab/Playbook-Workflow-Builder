import React from 'react'
import dynamic from 'next/dynamic'
import { connection_icon } from '@/icons'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function CAVATICAButton({ session_id }: { session_id?: string }) {
  const [disabled, setDisabled] = React.useState(false)
  return (
    <>
      <button
        className="bp4-button bp4-minimal"
      >
        <Icon
          title={session_id ? 'Disconnect from CAVATICA' : 'Connect to CAVATICA'}
          icon={connection_icon}
          className={disabled ? 'fill-gray-400' : session_id ? 'fill-red-600' : 'fill-black dark:fill-white'}
        />
      </button>
    </>
  )
}
