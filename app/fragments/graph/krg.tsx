import React from 'react'
import * as dict from '@/utils/dict'
import { MetaNode } from '@/spec/metanode'
import krg from '@/app/krg'
import useSWR from 'swr'
import dynamic from 'next/dynamic'

const UserIdentity = dynamic(() => import('@/app/fragments/graph/useridentity'))

export default function useKRG({ session_id }: { session_id?: string }) {
  const { data: suggestions, error } = useSWR(`/api/suggest`, { refreshInterval: 0 })
  const [krg_, setKrg_] = React.useState({ krg })
  React.useEffect(() => {
    if (!suggestions) return
    for (const suggestion of dict.values(suggestions)) {
      let OutputNode = krg.getDataNode(suggestion.output)
      if (OutputNode === undefined) {
        OutputNode = MetaNode(suggestion.output)
          .meta({
            label: `${suggestion.output} (Suggestion)`,
            description: `A data type, suggested as part of ${suggestion.name}`,
            pagerank: -100,
          })
          .codec<any>()
          .view((props) => {
            return <div>This data type was suggested as part of {suggestion.name}</div>
          })
          .build()
        krg.add(OutputNode)
        setKrg_({ krg })
      }
      let ProcessNode = krg.getProcessNode(suggestion.name)
      if (ProcessNode === undefined) {
        const ProcessNode = MetaNode(suggestion.name)
          .meta({
            label: `${suggestion.name} (Suggestion)`,
            description: suggestion.description,
            pagerank: -100,
          })
          .inputs(suggestion.inputs ?
              dict.init((suggestion.inputs as string).split(',').map((spec: string, ind: number) =>
                ({ key: ind.toString(), value: krg.getDataNode(spec) }))
                .filter(({ key, value }) => !!value))
              : {} as any)
          .output(OutputNode)
          .prompt((props) => {
            return <div>
              <p>{suggestion.description}</p>
              <p>This was suggested by {suggestion.user ? <UserIdentity user={suggestion.user} /> : <>a playbook workflow builder user</>}.</p>
            </div>
          })
          .story(props => ({ abstract: `It is suggested that "${suggestion.description}" be applied to the inputs: ${suggestion.inputs} to get a ${OutputNode.meta.label}.` }))
          .build()
        krg.add(ProcessNode)
        setKrg_({ krg })
      }
    }
  }, [suggestions])
  return krg_.krg
}
