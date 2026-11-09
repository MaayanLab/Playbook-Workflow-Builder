import React from 'react'
import { DataMetaNode, PromptMetaNode } from '@/spec/metanode'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { func_icon } from '@/icons'
import dynamic from 'next/dynamic'
import { Metapath, useMetapathInputs } from '@/app/fragments/metapath'
import type KRG from '@/core/KRG'
import { useStory } from '@/app/fragments/story'
import classNames from 'classnames'
import { AbstractPart, FigureCaption } from './story'
import { useTRPCClient } from '@/app/trpc/client'
import { useQueryClient } from '@tanstack/react-query'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function Prompt({ session_id, krg, processNode, outputNode, output, id, head, status, onRebase, collapseKey }: { session_id?: string, krg: KRG, processNode: PromptMetaNode, outputNode?: DataMetaNode, output: any, id: string, head: Metapath, status: string | undefined, onRebase: (headId: string) => void, collapseKey: string }) {
  const trpc = useTRPCClient()
  const queryClient = useQueryClient()
  const { data: inputs, error } = useMetapathInputs({ krg, head })
  const story = useStory()
  const [expanded, setExpanded] = React.useState(() => {
    if (typeof window === 'undefined') return true
    return window.sessionStorage.getItem(`pwb:report-collapse:${collapseKey}`) !== 'false'
  })
  React.useLayoutEffect(() => {
    window.sessionStorage.setItem(`pwb:report-collapse:${collapseKey}`, String(expanded))
  }, [collapseKey, expanded])
  const data = React.useMemo(() => {
    // invalidate data if it no longer matches the codec
    //  which could happen if an older verison of the metanode
    //  was used originally
    try {
      if (head.process.data !== null) {
        return processNode.codec.decode(head.process.data.value)
      }
    } catch (e) {}
  }, [head])
  const Component = processNode.prompt
  return (
    <div className={classNames('collapse collapse-arrow text-black dark:text-white', { 'collapse-open': expanded, 'bg-yellow-300 dark:bg-yellow-600': output === null, 'bg-green-400 dark:bg-green-700': output !== null })}>
      <div className="collapse-title flex flex-col gap-2 cursor-pointer" role="button" tabIndex={0} aria-expanded={expanded}
        onKeyDown={event => { if (event.target === event.currentTarget && (event.key === 'Enter' || event.key === ' ')) { event.preventDefault(); event.currentTarget.click() } }}
        onClick={event => { if (!(event.target as HTMLElement).closest('a')) setExpanded(value => !value) }}>
        <div className="flex flex-row gap-2 z-10">
          <Icon icon={processNode.meta.icon || func_icon} title={processNode.meta.label} className="fill-black dark:fill-white" />
          <h2 className="bp6-heading">{processNode.meta.label || processNode.spec}</h2>
        </div>
        <AbstractPart id={head.id} story={story} />
      </div>
      {error ? <div className="alert alert-error prose max-w-none">{error.toString()}</div> : null}
      <div className="collapse-content">
        {inputs !== undefined && array.intersection(dict.keys(processNode.inputs), dict.keys(inputs)).length === dict.keys(processNode.inputs).length ?
          <Component
            session_id={session_id}
            data={data}
            inputs={inputs}
            output={outputNode?.spec === processNode.output.spec ? output : undefined}
            submit={async (data) => {
              const res = await trpc.fpprg.rebase.mutate({
                fpl_id: id,
                process_id: head.process.id,
                process: {
                  type: head.process.type,
                  data: {
                    type: processNode.spec,
                    value: processNode.codec.encode(data),
                  },
                  inputs: head.process.inputs,
                },
              }, { context: { session_id } })
              // Fill the new query before changing the URL, keeping the current
              // report (and its scroll container) mounted throughout the edit.
              await queryClient.fetchQuery({
                queryKey: ['fpprg-fpl', session_id, res.head],
                queryFn: () => trpc.fpprg.fpl.query(res.head, { context: { session_id } }),
                staleTime: Infinity,
              })
              onRebase(res.head)
              const url = new URL(window.location.href)
              url.pathname = `${session_id ? `/session/${session_id}` : ''}/report/${res.head}`
              url.hash = ''
              window.history.pushState(null, '', url)
            }}
          />
          : <div className="prose max-w-none">Waiting for input</div>}
        {/* <FigureCaption id={head.id} story={story} /> */}
      </div>
      {status ? (
        <div className="alert shadow-lg place-content-start align-middle">
          <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" className="stroke-info shrink-0 w-6 h-6"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth="2" d="M13 16h-1v-4h-1m1-4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"></path></svg>
          <code className="prose max-w-none whitespace-pre-line">{status}</code>
        </div>
      ) : null}
      {outputNode && outputNode.spec === 'Error' && output ? outputNode.view(output) : null}
    </div>
  )
}
