import krg from '@/app/krg'
import db from '@/app/db'
import * as dict from '@/utils/dict'
import type { NextApiRequest, NextApiResponse } from 'next'
import dynamic from 'next/dynamic'
import { z } from 'zod'
import { MetaNode } from '@/spec/metanode'

const UserIdentity = dynamic(() => import('@/app/fragments/graph/useridentity'))

const BodyType = z.object({
  name: z.string(),
  inputs: z.string(),
  output: z.string(),
  user: z.string(),
  description: z.string(),
})

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (req.method === 'GET')  {
      res.status(200).json(await db.objects.suggestion.findMany())
    } else if (req.method === 'POST') {
      const suggestion = BodyType.parse(JSON.parse(req.body))
      // add the suggested KRG node(s)
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
      const ProcessNode = MetaNode(suggestion.name)
        .meta({
          label: `${suggestion.name} (Suggestion)`,
          description: suggestion.description,
          pagerank: -100,
        })
        .inputs(suggestion.inputs ?
            dict.init(suggestion.inputs.split(',').map((spec, ind) =>
              ({ key: ind.toString(), value: krg.getDataNode(spec) })))
            : {} as any)
        .output(OutputNode)
        .prompt((props) => {
          return <div>
            <p>{suggestion.description}</p>
            <p>This was suggested by {suggestion.user ? <UserIdentity user={suggestion.user} /> : <>a playbook partnership user</>}.</p>
          </div>
        })
        .story(props => `It is suggested that "${suggestion.description}" be applied to the inputs: ${suggestion.inputs} to get a ${OutputNode.meta.label}.`)
        .build()
      krg.add(ProcessNode)
      await db.objects.suggestion.create({ data: suggestion })
      res.status(200).end()
    } else {
      throw new Error('Unsupported method')
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
