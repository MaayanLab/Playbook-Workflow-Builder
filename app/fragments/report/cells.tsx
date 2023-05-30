import React from 'react'
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import type { Metapath } from '@/app/fragments/metapath'
import { StoryProvider } from '@/app/fragments/story'
import { useAPIMutation, useAPIQuery } from '@/core/api/client'
import { UserPlaybook, UpdateUserPlaybook, DeleteUserPlaybook, PublishUserPlaybook } from '@/app/api/client'
import * as dict from '@/utils/dict'
import { useRouter } from 'next/router'

const Introduction = dynamic(() => import('@/app/fragments/report/introduction'))
const Cell = dynamic(() => import('@/app/fragments/report/cell'))

export default function Cells({ krg, id }: { krg: KRG, id: string }) {
  const router = useRouter()
  const { data, error } = useAPIQuery(UserPlaybook, { id }, {
    keepPreviousData: true,
  })
  const { trigger: updateUserPlaybook } = useAPIMutation(UpdateUserPlaybook)
  const { trigger: publishUserPlaybook } = useAPIMutation(PublishUserPlaybook)
  const { trigger: deleteUserPlaybook } = useAPIMutation(DeleteUserPlaybook)
  const [cellMetadata, setCellMetadata] = React.useState({} as Record<string, Exclude<Metapath['cell_metadata'], null>>)
  const [playbookMetadata, setPlaybookMetadata] = React.useState({
    id: '',
    title: '',
    description: '',
    summary: 'auto',
    gpt_summary: '',
  } as Exclude<Metapath['playbook_metadata'], null>)
  const [userPlaybook, setUserPlaybook] = React.useState(undefined as undefined | { public: boolean })
  const [updateRequired, setUpdateRequired] = React.useState(false)
  React.useEffect(() => {
    if (!data) return
    setUpdateRequired(
      playbookMetadata?.id !== data.metapath[data.metapath.length-1].playbook_metadata?.id
      || data.metapath.some(cell => cell.cell_metadata?.id !== cellMetadata[cell.id]?.id)
    )
  }, [data, cellMetadata, playbookMetadata, userPlaybook])
  React.useEffect(() => {
    if (!data) return
    const {
      id = '',
      title = '',
      description = '',
      summary = 'auto',
      gpt_summary = '',
    } = data.metapath[data.metapath.length-1].playbook_metadata ?? {}
    setPlaybookMetadata({ id, title, description, summary, gpt_summary })
    setCellMetadata(dict.init(data.metapath.map((element, index) => {
      const {
        id = '',
        label = '',
        description = '',
        process_visible = false,
        data_visible = index+1 !== data.metapath.length-1,
      } = element.cell_metadata ?? {}
      return { key: element.id, value: { id, label, description, process_visible, data_visible } }
    })))
    setUserPlaybook(data.userPlaybook)
  }, [data])
  if (!data || !playbookMetadata) return null
  return (
    <div className="flex flex-col py-4 gap-2">
      <StoryProvider krg={krg} metapath={data.metapath}>
        <Introduction
          id={id}
          error={error}
          playbookMetadata={playbookMetadata}
          userPlaybook={userPlaybook}
          setPlaybookMetadata={setPlaybookMetadata}
          updateRequired={updateRequired}
          toggleSave={() => {
            if (updateRequired || !userPlaybook) {
              const { id: playbookMetadataId, ...playbook_metadata } = playbookMetadata
              const cell_metadata = dict.init(dict.items(cellMetadata).map(({ key, value }) => {
                const { id, ...meta } = value
                return { key, value: id ? { id } : meta }
              }))
              updateUserPlaybook({
                query: { id },
                body: {
                  user_playbook: { public: userPlaybook?.public || false },
                  playbook_metadata: playbookMetadataId ? { id: playbookMetadataId } : playbook_metadata,
                  cell_metadata,
                },
              }).then(id => {
                setUserPlaybook({ public: userPlaybook?.public || false })
                setPlaybookMetadata(metadata => ({ ...metadata, id: data.metapath[data.metapath.length-1].playbook_metadata?.id || '' }))
                setUpdateRequired(false)
                router.push(`/report/${id}`, undefined, { shallow: true, scroll: false })
              })
            } else {
              deleteUserPlaybook({
                query: { id },
                body: {},
              }).then(id => {
                setUserPlaybook(undefined)
                setUpdateRequired(false)
              })
            }
          }}
          togglePublic={() => {
            if (!updateRequired && userPlaybook) {
              const publicPlaybook = !userPlaybook.public
              publishUserPlaybook({
                query: { id },
                body: { public: publicPlaybook },
              }).then(id => {
                setUserPlaybook({ public: publicPlaybook })
                setUpdateRequired(false)
              })
            }
          }}
        />
        {(data.metapath||[]).filter(head => head.id in cellMetadata).map(head => (
          <Cell
            key={`${head.id}-${cellMetadata[head.id]?.process_visible}-${cellMetadata[head.id]?.data_visible}`}
            krg={krg}
            id={id}
            head={head}
            cellMetadata={cellMetadata}
            setCellMetadata={setCellMetadata}
          />
        ))}
      </StoryProvider>
    </div>
  )
}