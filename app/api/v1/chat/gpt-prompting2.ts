import { API } from '@/spec/api'
import { z } from 'zod'
import krg from '@/app/krg'
import * as dict from '@/utils/dict'
import type { ProcessMetaNode } from '@/spec/metanode'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { UnauthorizedError } from '@/spec/error'

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

async function openaiCompletions({ messages = [], model = 'gpt-3.5-turbo', temperature = 0.1, ...kwargs }: Partial<z.TypeOf<typeof OpenAICompletionsRequest>>) {
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

function strip(s: string) {
  return s.replace(/^\s+/g, '').replace(/\s+$/g, '')
}

type Component = { id: number, inputs?: Record<string, string>, data?: string, type: string, description: string }

function story(components: Component[]) {
  return components.map(component => {
    const componentProc = krg.getProcessNode(component.type)
    // todo include inputs
    return componentProc.story({})
  }).join(' ')
}

async function findAnswers(messages: { role: string, content: string }[]) {
  // find all previous components from message history
  const previousComponents = messages.reduce((components, { role, content }) => {
    if (role !== 'assistant') return components
    const m = /```component\n(.+)\n```/gm.exec(content)
    if (!m) return components
    return [...components, JSON.parse(m[1])]
  }, [] as Component[])
  const relevantMessages = messages.reduce((M, { role, content }) => {
    if (role === 'assistant') return M
    const m = /```(component|suggestions)\n(.+)\n```/gm.exec(content)
    if (m) return M
    return [...M, { role, content }]
  }, [] as {role: string, content: string}[])
  const ctx = { maxId: Math.max(0, ...previousComponents.map(component => component.id)) }

  const results: { role: 'assistant', content: string }[] = []

  for (let i = 0; i < 3; i++) {
    // construct possible next components
    const nextComponents = dict.init([
      ...krg.getNextProcess('')
        .map(proc => ({ id: ++ctx.maxId, type: proc.spec, description: proc.meta.description })),
      ...previousComponents.flatMap(component => {
        const componentProc = krg.getProcessNode(component.type)
        return krg.getNextProcess(componentProc.output.spec)
          .filter(proc => proc.spec !== componentProc.spec)
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
        ...relevantMessages,
        {
          role: 'user',
          content: [
            ...(previousComponents.length > 0 ? [
              'The current workflow description is:\n',
              '```',
              story(previousComponents),
              '```',
            ] : []),
            'Considering my messages up to this point, and the following potential next workflow components:',
            '```json',
            JSON.stringify(Object.values(nextComponents).map(({ type: _, ...component }) => component)),
            '```',
            'Respond with only a JSON serialized document which satisfies the typescript type:',
            '```ts',
            `type Step = {`,
            ` /* an id corresponding to a "next" component above */ "id": number,`,
            ` /* part of the user's messages relevant to the component */ "term": string }`,
            `const Response: {`,
            ` /* response to the user with your thoughts */ "response": string,`,
            ` /* at most the best 3 good choices for the possible next step to accomplish the user's instructions */ "choices": Step[] } =`,
          ].join('\n'),
        },
      ]
    })
    const content = strip(res.choices[0].message.content)
    try {
      const m = /\{.+\}/gms.exec(content)
      const contentParsed = z.object({
        response: z.string(),
        choices: z.union([z.object({ id: z.number(), term: z.string() }), z.array(z.object({ id: z.number(), term: z.string() }))]).transform(choice => Array.isArray(choice) ? choice : [choice]),
      }).parse(JSON.parse(m ? m[0] : ''))
      const possibleComponents = contentParsed.choices.map(choice => {
        const component = nextComponents[choice.id]
        if (component.type.includes('Input')) {
          Object.assign(component, {data: choice.term})
        }
        return component
      })
      if (contentParsed.response) {
        results.push({ role: 'assistant', content: contentParsed.response })
      }
      if (possibleComponents.length === 1) {
        const [possibleComponent] = possibleComponents
        previousComponents.push(possibleComponent)
        results.push({ role: 'assistant', content: '```component\n' + JSON.stringify(possibleComponent) + '\n```' })
      } else if (possibleComponents.length > 1) {
        results.push({ role: 'assistant', content: '```suggestions\n' + JSON.stringify(possibleComponents) + '\n```' })
        break
      } else {
        break
      }
    } catch (e) {
      results.push({ role: 'assistant', content: `I'm not sure, please try rephrasing your question.` })
      break
    }
  }
  return results
}

export const ChatGPTPrompting2 = API('/api/v1/chat/gpt-prompting2')
  .query(z.object({}))
  .body(z.object({
    messages: z.array(z.object({
      role: z.string(),
      content: z.string(),
    })),
  }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }

    const newMessages: typeof inputs.body.messages = []
    try {
      newMessages.push(...(await findAnswers(inputs.body.messages)))
    } catch (e) {
      newMessages.push({ role: 'assistant', content: `Error: ${(e as any).toString()}` })
    }
    return newMessages
  })
  .build()
