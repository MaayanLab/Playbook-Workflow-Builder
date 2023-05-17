import React from 'react'
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import type { Metapath } from '@/app/fragments/metapath'
import { useSWRImmutableSticky } from '@/utils/use-sticky'
import { StoryProvider } from '@/app/fragments/story'
import useAsyncEffect from 'use-async-effect'
import * as dict from '@/utils/dict'

const Introduction = dynamic(() => import('@/app/fragments/report/introduction'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))

export type CellMetadata = {
  label: string, description: string,
  gpt_summary?: string, summary: 'auto' | 'gpt' | 'manual',
  processVisible: boolean, dataVisible: boolean,
}

export type Playbook = {
  id?: string
  public?: string,
  update_required: boolean,
}

export default function Cells({ krg, id }: { krg: KRG, id: string }) {
  const { data: metapath, error } = useSWRImmutableSticky<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
  const [metadata, setMetadata_] = React.useState({} as Record<string, CellMetadata>)
  const [playbook, setPlaybook] = React.useState({
    update_required: false,
  } as Playbook)
  const setMetadata = React.useCallback(
    (cb: (currentValue: Record<string, CellMetadata>) => Record<string, CellMetadata>) => {
      setMetadata_(cb)
      setPlaybook(p => ({ ...p, update_required: p.id !== undefined }))
  }, [setMetadata_]) as React.Dispatch<React.SetStateAction<Record<string, CellMetadata>>>
  React.useEffect(() => {
    if (!metapath) return
    setMetadata(() => dict.init(metapath.map((element, index) => {
      const {
        label = '',
        description = '',
        processVisible = false,
        dataVisible = index+1 !== metapath.length-1,
      } = element.metadata ?? {}
      const summary = 'auto'
      return { key: element.id, value: { label, description, summary, processVisible, dataVisible } }
    })))
  }, [metapath])
  if (!metapath || !metadata[id]) return null
  return (
    <div className="flex flex-col py-4 gap-2">
      <StoryProvider krg={krg} metapath={metapath}>
        <Introduction
          id={id}
          error={error}
          metadata={metadata}
          setMetadata={setMetadata}
          playbook={playbook}
          setPlaybook={setPlaybook}
        />
        {(metapath||[]).filter(head => head.id in metadata).map(head => (
          <Cell
            key={`${head.id}-${metadata[head.id].processVisible}-${metadata[head.id].dataVisible}`}
            krg={krg}
            id={id}
            head={head}
            metadata={metadata}
            setMetadata={setMetadata}
          />
        ))}
      </StoryProvider>
    </div>
  )
}