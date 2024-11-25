import React from 'react'
import dynamic from 'next/dynamic'
import { export_icon } from '@/icons'

const Icon = dynamic(() => import('@/app/components/icon'))

const Bp5Popover = dynamic(() => import('@blueprintjs/core').then(({ Popover }) => Popover))
const Bp5Menu = dynamic(() => import('@blueprintjs/core').then(({ Menu }) => Menu))
const Bp5MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))

export default function ExportButton({ session_id, id, metadata }: { session_id?: string, id: string, metadata: { title: string, description: string | undefined } }) {
  return (
    <Bp5Popover
      className={'cursor-pointer'}
      content={
        <Bp5Menu>
          <Bp5MenuItem
            icon="document"
            text="Playbook JSON"
            href={`${session_id ? `/api/socket/${session_id}` : ''}/api/db/fpl/${id}/export`}
            download={`playbook-${id}.json`}
          />
          <Bp5MenuItem
            icon="document"
            text="Playbook with Output JSON"
            href={`${session_id ? `/api/socket/${session_id}` : ''}/api/db/fpl/${id}/output/export`}
            download={`playbook-${id}.json`}
          />
          <Bp5MenuItem
            icon="document"
            text="BioCompute Object (BCO) JSON"
            href={`${session_id ? `/api/socket/${session_id}` : ''}/api/bco/${id}?metadata=${encodeURIComponent(JSON.stringify(metadata))}`}
            download={`bco-${id}.json`}
          />
          <Bp5MenuItem
            icon="document"
            text="Common Workflow Language (CWL) Bundle"
            href={`${session_id ? `/api/socket/${session_id}` : ''}/api/v1/cwl/${id}`}
            download={`cwl-${id}.zip`}
          />
          <Bp5MenuItem
            icon="document"
            text="Research Object Crate (RO-Crate)"
            href={`${session_id ? `/api/socket/${session_id}` : ''}/api/v1/ro-crate/${id}`}
            download={`ro-crate-metadata.json`}
          />
        </Bp5Menu>
      }
      placement="bottom"
    >
      <Icon
        icon={export_icon}
        className={'fill-black dark:fill-white'}
        title="Export Workflow"
      />
    </Bp5Popover>
  )
}
