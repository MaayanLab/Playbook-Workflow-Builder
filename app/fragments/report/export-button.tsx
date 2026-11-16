import React from 'react'
import dynamic from 'next/dynamic'
import { export_icon } from '@/icons'
import usePublicUrl from '@/utils/next-public-url'
import { useExport } from '@/app/fragments/report/use-export'
import { useSessionWithId } from '@/app/extensions/auth/hooks'
import * as Auth from '@/app/extensions/auth-client'

const Icon = dynamic(() => import('@/app/components/icon'))

const BpPopoverNext = dynamic(() => import('@blueprintjs/core').then(({ PopoverNext }) => PopoverNext))
const Bp5Menu = dynamic(() => import('@blueprintjs/core').then(({ Menu }) => Menu))
const Bp5MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))

export default function ExportButton({ session_id, id, metadata }: { session_id?: string, id: string, metadata: { title: string, description: string | undefined } }) {
  const publicUrl = usePublicUrl({ absolute: true })
  // these can take a while to produce, so they're produced in the background (see use-export)
  //  & kept as the user's uploads, which requires signing in
  const startExport = useExport(session_id)
  const userSession = useSessionWithId()
  const signedIn = !!userSession.data?.user
  const signInTip = 'Sign in to produce this export, it will be kept with your uploads'
  return (
    <BpPopoverNext
      className={'cursor-pointer'}
      content={
        <Bp5Menu>
          {!signedIn && userSession.status !== 'loading' ?
            <Bp5MenuItem
              icon="log-in"
              text="Sign in to export as CWL, RO-Crate, LaTeX or PDF"
              onClick={() => { Auth.signIn() }}
            />
            : null}
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
          <div className="tooltip block text-left" data-tip={signedIn ? "Download a bundle of CWL documents which can be used to execute this workflow with a CWL runner" : signInTip}>
            <Bp5MenuItem
              disabled={!signedIn}
              icon="document"
              text="Common Workflow Language (CWL) Bundle"
              onClick={() => startExport({ fpl_id: id, format: 'cwl-zip', label: 'CWL Bundle' })}
            />
          </div>
          <div className="tooltip block text-left" data-tip={signedIn ? "Download RO-Crate object bundle for this workflow" : signInTip}>
            <Bp5MenuItem
              disabled={!signedIn}
              icon="document"
              text="Research Object Crate (RO-Crate)"
              onClick={() => startExport({ fpl_id: id, format: 'ro-crate-zip', label: 'Research Object Crate' })}
            />
          </div>
          <div className="tooltip block text-left" data-tip="Download RO-Crate JSON object metadata for this workflow">
            <Bp5MenuItem
              icon="document"
              text="Research Object Crate (RO-Crate) (Detached)"
              href={`${session_id ? `/api/socket/${session_id}` : ''}/api/v1/ro-crate/${id}?format=json`}
              download={`ro-crate-metadata.json`}
            />
          </div>
          <div className="tooltip block text-left" data-tip={signedIn ? "Download a latex rendering of the playbook for creating a paper. This feature is in early BETA and is currently incomplete." : signInTip}>
            <Bp5MenuItem
              disabled={!signedIn}
              icon="document"
              text="LaTeX Bundle (BETA)"
              onClick={() => startExport({ fpl_id: id, format: 'tex-zip', metadata, label: 'LaTeX Bundle' })}
            />
          </div>
          {publicUrl ? <div className="tooltip block text-left" data-tip={signedIn ? "Download a latex rendering of the playbook for creating a paper. This feature is in early BETA and is currently incomplete." : signInTip}>
            <Bp5MenuItem
              disabled={!signedIn}
              icon="document"
              text="LaTeX Bundle in Overleaf (BETA)"
              onClick={() => startExport({ fpl_id: id, format: 'tex-zip', metadata, label: 'LaTeX Bundle', then: 'overleaf' })}
            />
          </div>
          : null}
          <div>
            <div className="tooltip block text-left" data-tip={signedIn ? "Download a PDF rendering of the playbook report. This feature is in early BETA and is currently incomplete." : signInTip}>
              <Bp5MenuItem
                disabled={!signedIn}
                icon="document"
                text="Report PDF (BETA)"
                onClick={() => startExport({ fpl_id: id, format: 'tex-pdf', metadata, label: 'Report PDF' })}
              />
            </div>
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
    </BpPopoverNext>
  )
}
