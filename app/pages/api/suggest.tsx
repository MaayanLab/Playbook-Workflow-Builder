import krg from '@/app/krg'
import db from '@/app/kvdb'
import * as dict from '@/utils/dict'
import type { NextApiRequest, NextApiResponse } from 'next'
import { z } from 'zod'
import { MetaNode } from '@/spec/metanode'
import { v4 as uuidv4 } from 'uuid'
import sub from 'subleveldown'

const dbSuggest = sub(db, 'suggest')

const BodyType = z.object({
  name: z.string(),
  inputs: z.string(),
  output: z.string(),
  author_name: z.string(),
  author_email: z.string(),
  author_org: z.string(),
  description: z.string(),
})

export default async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    if (req.method === 'GET')  {
      const kvdb: Record<string, string> = {}
      for await (const [key, value] of dbSuggest.iterator() as any) {
        kvdb[key] = value.toString()
      }
      res.status(200).json(kvdb)
    } else if (req.method === 'POST') {
      const suggestion = BodyType.parse(JSON.parse(req.body))
      // add the suggested KRG node(s)
      let OutputNode = krg.getDataNode(suggestion.output)
      if (OutputNode === undefined) {
        OutputNode = MetaNode.createData(suggestion.output)
          .meta({
            label: suggestion.output,
            description: `A data type, suggested as part of ${suggestion.name}`,
          })
          .codec<any>()
          .view((props) => {
            return <div>This data type was suggested as part of {suggestion.name}</div>
          })
          .build()
        krg.add(OutputNode)
      }
      const ProcessNode = MetaNode.createProcess(suggestion.name)
        .meta({
          label: suggestion.name,
          description: suggestion.description,
        })
        .inputs(suggestion.inputs ?
            dict.init(suggestion.inputs.split(',').map((spec, ind) =>
              ({ key: ind.toString(), value: krg.getDataNode(spec) })))
            : {} as any)
        .output(OutputNode)
        .prompt((props) => {
          return <div>This was suggested by {suggestion.author_name} &lt;{suggestion.author_email}&gt; ({suggestion.author_org})</div>
        })
        .build()
      krg.add(ProcessNode)
      // register the suggestion
      await dbSuggest.put(`${(new Date()).toString()}-${uuidv4()}`, JSON.stringify(suggestion))
      res.status(200).end()
    } else {
      throw new Error('Unsupported method')
    }
  } catch (e) {
    console.error(e)
    res.status(500).end((e as Error).toString())
  }
}
