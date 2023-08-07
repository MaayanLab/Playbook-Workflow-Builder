import React from 'react'
import dynamic from 'next/dynamic'
import { download_icon } from '@/icons'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function DownloadButton({ id }: { id: string }) {
  return (
    <a
      className="bp4-button bp4-minimal"
      href={`/api/db/fpl/${id}/output`}
      download={`${id}.json`}
    >
      <Icon
        icon={download_icon}
        className={'fill-black dark:fill-white'}
        title="Download Report"
      />
    </a>
  )
}
