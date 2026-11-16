/**
 * User suggestions (see the `suggestion` table) are represented in the KRG as
 *  placeholder nodes, so that playbooks containing suggested steps can be shown,
 *  saved and exported. Shared by the server (@/app/server-krg) and the browser
 *  (@/app/fragments/metapath).
 */
import dynamic from 'next/dynamic'
import type KRG from '@/core/KRG'
import type { Suggestion } from '@/db/schema'
import * as dict from '@/utils/dict'
import { MetaNode } from '@/spec/metanode'

const UserIdentity = dynamic(() => import('@/app/fragments/graph/useridentity'))

export type SuggestionSpec = Pick<Suggestion, 'name' | 'inputs' | 'output' | 'description' | 'user'>

/** Suggestions belong to the database behind a local or tunneled session. */
export const suggestionsQueryKey = (session_id?: string) => ['catalog-suggestions', session_id] as const

/**
 * A KRG including the given suggestions (`krg` itself is not modified)
 */
export function withSuggestions(krg: KRG, suggestions: SuggestionSpec[]): KRG {
  if (suggestions.length === 0) return krg
  const extended = krg.clone()
  for (const suggestion of suggestions) addSuggestion(extended, suggestion)
  return extended
}

/**
 * Add a suggestion's process (& output data type, unless it exists) to `krg`
 */
export function addSuggestion(krg: KRG, suggestion: SuggestionSpec) {
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
  }
  if (krg.getProcessNode(suggestion.name) === undefined) {
    const Output = OutputNode
    const ProcessNode = MetaNode(suggestion.name)
      .meta({
        label: `${suggestion.name} (Suggestion)`,
        description: suggestion.description,
        pagerank: -100,
      })
      .inputs(suggestion.inputs ?
          dict.init(suggestion.inputs.split(',').map((spec: string, ind: number) =>
            ({ key: ind.toString(), value: krg.getDataNode(spec) }))
            .filter(({ key, value }) => !!value))
          : {} as any)
      .output(Output)
      .prompt((props) => {
        return <div>
          <p>{suggestion.description}</p>
          <p>This was suggested by {suggestion.user ? <UserIdentity user={suggestion.user} /> : <>a playbook workflow builder user</>}.</p>
        </div>
      })
      .story(props => ({ abstract: `It is suggested that "${suggestion.description}" be applied to the inputs: ${suggestion.inputs} to get a ${Output.meta.label}.` }))
      .build()
    krg.add(ProcessNode)
  }
}
