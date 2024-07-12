import React from 'react'
import dynamic from 'next/dynamic'
import { import_icon } from '@/icons'
import { useRouter } from 'next/router'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function ImportButton({ session_id }: { session_id?: string }) {
  const router = useRouter()
  // TODO: error dialog on failure
  return (
    <>
      <input
        id="import-button"
        type="file"
        className="hidden"
        accept="application/json"
        onChange={async (evt) => {
          evt.preventDefault()
          const body = await new Promise<string | null>((resolve, reject) => {
            if (evt.target.files === null) return null
            const fileReader = new FileReader()
            fileReader.readAsText(evt.target.files[0], 'UTF-8')
            fileReader.addEventListener('load', () => resolve(fileReader.result as string))
            fileReader.addEventListener('error', (e) => reject(e))
          })
          const req = await fetch(`${session_id ? `/api/socket/${session_id}` : ''}/api/db/fpl`, {
            headers: {
              'Content-Type': 'application/json',
            },
            method: 'POST',
            body,
          })
          const res = await req.json()
          if (req.ok) {
            router.push(`${session_id ? `/session/${session_id}` : ''}/report/${res}`)
          } else {
            console.error(res)
          }
        }}
      />
      <label htmlFor="import-button" className="bp5-button bp5-minimal">
        <Icon
          icon={import_icon}
          className="fill-black dark:fill-white"
          title="Import workflow from Playbook JSON"
        />
      </label>
    </>
  )
}
