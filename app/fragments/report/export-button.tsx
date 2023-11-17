import React from 'react'
import dynamic from 'next/dynamic'
import { export_icon } from '@/icons'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function ExportButton({ session_id, id }: { session_id?: string, id: string }) {
  return (
    <a
      className="bp5-button bp5-minimal"
      href={`${session_id ? `/api/socket/${session_id}` : ''}/api/db/fpl/${id}/output/export`}
      download={`${id}.json`}
    >
      <Icon
        icon={export_icon}
        className={'fill-black dark:fill-white'}
        title="Export Workflow"
      />
    </a>
  )
}
