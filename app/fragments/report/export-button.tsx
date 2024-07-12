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
          <div className="tooltip block text-left" data-tip="Download the workflow steps encoded as JSON">
            <Bp5MenuItem
              icon="document"
              text="Playbook JSON"
              href={`${session_id ? `/api/socket/${session_id}` : ''}/api/db/fpl/${id}/export`}
              download={`playbook-${id}.json`}
            />
          </div>
          <div className="tooltip block text-left" data-tip="Download the workflow steps along with the output encoded as JSON">
            <Bp5MenuItem
              icon="document"
              text="Playbook with Output JSON"
              href={`${session_id ? `/api/socket/${session_id}` : ''}/api/db/fpl/${id}/output/export`}
              download={`playbook-${id}.json`}
            />
          </div>
          <div className="tooltip block text-left" data-tip="Download the workflow's BCO JSON document compatible with IEEE 2791-2022">
              <Bp5MenuItem
              icon="document"
              text="BioCompute Object (BCO) JSON"
              href={`${session_id ? `/api/socket/${session_id}` : ''}/api/bco/${id}?metadata=${encodeURIComponent(JSON.stringify(metadata))}`}
              download={`bco-${id}.json`}
            />
          </div>
          <div className="tooltip block text-left" data-tip="Download a bundle of CWL documents which can be used to execute this workflow with a CWL runner">
            <Bp5MenuItem
              icon="document"
              text="Common Workflow Language (CWL) Bundle"
              href={`${session_id ? `/api/socket/${session_id}` : ''}/api/v1/cwl/${id}`}
              download={`cwl-${id}.zip`}
            />
          </div>
          </Bp5Menu>
      }
      placement="bottom"
    >
      <Icon
        icon={export_icon}
        className={'fill-black dark:fill-white'}
        title="Export workflow to several formats (click for options)"
      />
    </Bp5Popover>
  )
}
