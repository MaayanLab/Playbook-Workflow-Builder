import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import type { Metapath } from '@/app/fragments/metapath'
import { useSWRImmutableSticky } from '@/utils/use-sticky'
import { StoryProvider } from '@/app/fragments/report/story'

const Introduction = dynamic(() => import('@/app/fragments/report/introduction'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))

export default function Cells({ krg, id }: { krg: KRG, id: string }) {
  const { data: metapath, error } = useSWRImmutableSticky<Array<Metapath>>(id ? `/api/db/fpl/${id}` : undefined)
  if (!metapath) return null
  return (
    <div className="flex flex-col py-4 gap-2">
      <StoryProvider krg={krg} metapath={metapath}>
        <Introduction id={id} error={error} />
        {(metapath||[]).map((head, index) => (
          <Cell key={head.id} krg={krg} id={id} head={head} defaultCollapse={index+1 !== metapath?.length} />
        ))}
      </StoryProvider>
    </div>
  )
}