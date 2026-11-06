/**
 * Register user suggestions (see the `suggestion` table) as placeholder nodes in the
 *  server's KRG, server-side code (e.g. playbook saving & exports) looks up the
 *  process & data types of every step in a playbook, including suggested ones.
 */
import dynamic from 'next/dynamic'
import db from '@/app/db'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'
import { MetaNode } from '@/spec/metanode'

const UserIdentity = dynamic(() => import('@/app/fragments/graph/useridentity'))

export async function registerSuggestions() {
  const suggestions = await db.selectFrom('suggestion').selectAll().execute()
  for (const key in suggestions) {
    const suggestion = suggestions[key]
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
    let ProcessNode = krg.getProcessNode(suggestion.name)
    if (ProcessNode === undefined) {
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
    }
  }
}
