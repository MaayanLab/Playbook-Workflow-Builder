import React from 'react'
import dynamic from 'next/dynamic'
import { import_icon } from '@/icons'
import { useExRouter } from '@/app/fragments/ex-router'
import { useTRPCClient } from '@/app/trpc/client'
import { z } from 'zod'
import { DataC, IdOrPlaybookMetadataC, IdOrProcessC } from '@/core/FPPRG'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function ImportButton({ session_id }: { session_id?: string }) {
  const router = useExRouter()
  const trpc = useTRPCClient()
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
          try {
            const input = z.object({
              data: z.record(z.string(), DataC).optional(),
              workflow: z.array(IdOrProcessC),
              metadata: IdOrPlaybookMetadataC.optional(),
            }).parse(JSON.parse(body ?? ''))
            const res = await trpc.fpprg.import.mutate(input, { context: { session_id } })
            router.push(`${session_id ? `/session/${session_id}` : ''}/report/${res}`, undefined, { shallow: true })
          } catch (error) {
            console.error(error)
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
