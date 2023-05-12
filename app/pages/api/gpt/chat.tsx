import krg from '@/app/krg'
import handler from '@/utils/next-rest'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import type { ProcessMetaNode } from '@/spec/metanode'

function strip(s: string) {
  return s.replace(/^\s+/g, '').replace(/\s+$/g, '')
}

const OpenAICompletionsRequest = z.object({
  model: z.string(),
  temperature: z.number(),
  stop: z.array(z.string()),
  messages: z.array(z.object({
    role: z.string(),
    content: z.string(),
  })),
})

const OpenAICompletionsResponse = z.object({
  choices: z.array(z.object({
    message: z.object({
      role: z.string(),
      content: z.string(),
    })
  }))
})

async function openaiCompletions({ messages = [], model = 'gpt-3.5-turbo', temperature = 0.0, ...kwargs }: Partial<z.TypeOf<typeof OpenAICompletionsRequest>>) {
  console.debug(`send ${JSON.stringify(messages)}`)
  const req = await fetch(`https://api.openai.com/v1/chat/completions`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${process.env.OPENAI_API_KEY}`,
    },
    body: JSON.stringify({ messages, model, temperature, ...kwargs }),
  })
  const res = await req.json()
  console.debug(`recv ${JSON.stringify(res)}`)
  return OpenAICompletionsResponse.parse(res)
}

function nodeDesc(node: ProcessMetaNode, { inputs, output }: { inputs?: Record<string, unknown>, output?: string }) {
  if (inputs !== undefined && node.story !== undefined) {
    try {
      return `${node.meta.label}: ${node.meta.description}. ${node.story({
        inputs: inputs !== undefined ? inputs : {},
        output,
      })}`
    } catch (e) {}
  }
  return `${node.meta.label}: ${node.meta.description}`
}

type Component = { id: number, inputs?: Record<string, string>, data?: string, type: string, description: string }

async function findAnswers(messages: { role: "user" | "system" | "assistant", content: string }[]) {
  // find all previous components from message history
  const previousComponents = messages.reduce((components, { role, content }) => {
    if (role !== 'assistant') return components
    const m = /```component\n(.+)\n```/gm.exec(content)
    if (!m) return components
    return [...components, JSON.parse(m[1])]
  }, [] as Component[])
  // construct possible next components
  const ctx = { maxId: Math.max(0, ...previousComponents.map(component => component.id)) }
  const nextComponents = dict.init([
    ...krg.getNextProcess('')
      .map(proc => ({ id: ++ctx.maxId, type: proc.spec, description: proc.meta.description })),
    ...previousComponents.flatMap(component => {
      const componentProc = krg.getProcessNode(component.type)
      return krg.getNextProcess(componentProc.output.spec)
        .map(proc => {
          console.log(`${ctx.maxId} ${componentProc.spec} => ${componentProc.output.spec} => ${proc.spec}`)
          // todo use the same approach here as elsewhere to capture multi-inputs properly
          const inputs = {} as Record<string, unknown>
          const inputValues = {} as Record<string, unknown>
          dict.items(proc.inputs).forEach(({ key, value }) => {
            if (Array.isArray(value)) {
              inputs[key] = [component.id]
              inputValues[key] = [component.data]
            }
            else {
              inputs[key] = component.id
              inputValues[key] = component.data
            }
          })
          return { id: ++ctx.maxId, inputs, type: proc.spec, description: nodeDesc(proc, { inputs: inputValues }) }
        })
    })
  ].map(value => ({ key: value.id, value })))
  const res = await openaiCompletions({
    messages: [
      ...messages.reduce((M, { role, content }) => {
        if (role !== 'assistant') {
          return [...M, { role, content }]
        }
        const m = /```component\n(.+)\n```/gm.exec(content)
        if (!m) return M
        const component = JSON.parse(m[1])
        return [
          ...M,
          { role: 'assistant', content: JSON.stringify({ id: component.id, term: component.data }) },
          { role: 'assistant', content },
        ]
      }, [] as {role: string, content: string}[]),
      {
        role: 'user',
        content: [
          'Considering our conversation up to this point, especially my most recent message and the following potential next workflow components:',
          '```json',
          JSON.stringify(Object.values(nextComponents)),
          '```',
          'Respond with a JSON serialization which satisfies the typescript type:',
          '```ts',
          `const response: {`,
          ` /* an id corresponding to the component above which best satisfies the user's intent */ "id": number,`,
          ` /* part of the user's messages relevant to the component */ "term": string } | {`,
          ` /* only if there are no good choices */ "question": string } =`,
        ].join('\n'),
      },
    ]
  })
  const content = strip(res.choices[0].message.content)
  try {
    const m = /\{.+\}/gms.exec(content)
    const contentParsed = JSON.parse(m ? m[0] : '')
    if (typeof contentParsed.id === 'number') {
      if (contentParsed.id in nextComponents) {
        const component = nextComponents[contentParsed.id]
        if (component.type.includes('Input') && 'term' in contentParsed) Object.assign(component, {data: contentParsed.term})
        return '```component\n' + JSON.stringify(component) + '\n```'
      } else {
        return content
      }
    }
    return contentParsed.question
  } catch (e) {}
  return content
}

const BodyType = z.object({
  messages: z.array(z.object({
    role: z.enum(['system', 'user', 'assistant']),
    content: z.string(),
  })),
})

export default handler(async (req, res) => {
  const session = await getServerSessionWithId(req, res)
  if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
  if (req.method === 'HEAD') {
    res.status(200).end()
    return
  }
  if (req.method !== 'POST') throw new UnsupportedMethodError()
  const body = BodyType.parse(JSON.parse(req.body))
  const lastMessage = body.messages[body.messages.length-1]
  if (lastMessage.role !== 'user') throw new Error('Unexpected user input')
  const content = await findAnswers(body.messages)
  res.status(200).json({ role: 'assistant', content })
})
