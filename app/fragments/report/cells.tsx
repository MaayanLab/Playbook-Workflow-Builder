import React from 'react'
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import type { Metapath } from '@/app/fragments/metapath'
import { useSWRImmutableSticky } from '@/utils/use-sticky'
import { StoryProvider } from '@/app/fragments/story'
import useAsyncEffect from 'use-async-effect'

const Introduction = dynamic(() => import('@/app/fragments/report/introduction'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))

export type ReportMetadata = {
  title: string,
  description?: string,
  gpt_summary?: string,
  summary: 'auto' | 'manual' | 'gpt',
  collapsed: Record<number, boolean>,
  public: boolean,
}

export default function Cells({ krg, id }: { krg: KRG, id: string }) {
  const { data: metapath, error } = useSWRImmutableSticky<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
  const [metadata, setMetadata] = React.useState<ReportMetadata>({
    title: 'Playbook',
    summary: 'auto',
    collapsed: {},
    public: false,
  })
  useAsyncEffect(async (isMounted) => {
    setMetadata(({ gpt_summary: _, summary, ...metadata }) => ({
      ...metadata,
      summary: summary === 'gpt' ? 'auto' : summary,
    }))
    if (!metapath) return
    const published = (await import('@/app/public/playbooksDemo')).default
    if (!isMounted()) return
    published.filter((playbook) => id === playbook.id).forEach(playbook => {
      setMetadata((metadata) => ({
        ...metadata,
        summary: 'manual',
        title: playbook.label,
        collapsed: playbook.collapsed,
        description: playbook.description,
        gpt_summary: playbook.gpt_summary,
        public: true,
      }))
    })
  }, [id, metapath])
  if (!metapath) return null
  return (
    <div className="flex flex-col py-4 gap-2">
      <StoryProvider krg={krg} metapath={metapath}>
        <Introduction
          id={id}
          error={error}
          defaultMetadata={metadata}
        />
        {(metapath||[]).map((head, index) => (
          <Cell
            key={`${head.id}-${metadata.collapsed[index]}`}
            krg={krg}
            id={id}
            head={head}
            defaultCollapse={metadata.collapsed[index] === false ? false : index+1 !== metapath?.length}
          />
        ))}
      </StoryProvider>
    </div>
  )
}