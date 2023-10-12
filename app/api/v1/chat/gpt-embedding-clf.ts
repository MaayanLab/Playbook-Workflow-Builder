import { API } from '@/spec/api'
import krg from '@/app/krg'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError } from '@/spec/error'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import type { ProcessMetaNode } from '@/spec/metanode'
import gptEmbeddingClfDef from '@/app/public/chat/gpt-embedding-clf.json'
import { mean } from '@/utils/math'

async function gptEmbeddingClf({ embedding }: { embedding: number[] }) {
  const req = await fetch(`${gptEmbeddingClfDef.url}:predict`, {
    method: 'POST',
    body: JSON.stringify({
      instances: [embedding],
    }),
  })
  const res = await req.json()
  const { predictions: [predictions] } = z.object({ predictions: z.array(z.array(z.number())) }).parse(res)
  const scores: Record<string, number> = {}
  gptEmbeddingClfDef.vocab.forEach((workflow, i) => {
    scores[workflow] = predictions[i]
  })
  return scores
}

const OpenAIEmbeddingsRequest = z.object({
  model: z.string(),
  input: z.array(z.string()),
})

const OpenAIEmbeddingsResponse = z.object({
  data: z.array(z.object({
    embedding: z.array(z.number()),
  }))
})

async function openaiEmbedding({ input, model = 'text-embedding-ada-002', ...kwargs }: Partial<z.TypeOf<typeof OpenAIEmbeddingsRequest>>) {
  console.debug(`send ${JSON.stringify(input)}`)
  const req = await fetch(`https://api.openai.com/v1/embeddings`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${process.env.OPENAI_API_KEY}`,
    },
    body: JSON.stringify({ input, model, ...kwargs }),
  })
  const res = await req.json()
  // console.debug(`recv ${JSON.stringify(res)}`)
  return OpenAIEmbeddingsResponse.parse(res)
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

type Component = { id: number, inputs?: Record<string, number | number[]>, data?: string, type: string, description: string }

async function findAnswers(messages: { role: string, content: string }[]) {
  const lastMessage = array.findLast(messages, ({ role }) => role === 'user')
  if (!lastMessage) throw new Error('No user message')
  // obtain likelihood scores
  const [{ embedding }] = (await openaiEmbedding({ input: [lastMessage.content] })).data
  const scores = await gptEmbeddingClf({ embedding })
  // find all previous components from message history
  const previousComponents = messages.reduce((components, { role, content }) => {
    if (role !== 'assistant') return components
    const m = /```component\n(.+)\n```/gm.exec(content)
    if (!m) return components
    return [...components, JSON.parse(m[1])]
  }, [] as Component[])
  const ctx = { maxId: Math.max(0, ...previousComponents.map(component => component.id)) }
  const results: { role: 'assistant', content: string }[] = []
  // extend while high likelihood, otherwise offer suggestions
  for (let i = 0; i < 10; i++) {
    const nextComponents: {
      id: number,
      type: string,
      inputs?: Record<string, number | number[]>,
      description: string,
      likelihood: number,
    }[] = []
    if (previousComponents.length === 0) {
      krg.getNextProcess().forEach(proc => nextComponents.push({
        id: ++ctx.maxId,
        type: proc.spec,
        description: proc.meta.description,
        likelihood: scores[proc.spec] || 0,
      }))
    } else {
      const lastComponent = previousComponents[previousComponents.length-1]
      const lastComponentProc = krg.getProcessNode(lastComponent.type)
      krg.getNextProcess(lastComponentProc.output.spec)
        .forEach(proc => {
          // todo use the same approach here as elsewhere to capture multi-inputs properly
          const inputs = {} as Record<string, number | number[]>
          const inputValues = {} as Record<string, unknown>
          const likelihoods = [] as number[]
          dict.items(proc.inputs).forEach(({ key, value }) => {
            if (Array.isArray(value)) {
              inputs[key] = [lastComponent.id]
              inputValues[key] = [lastComponent.data]
            }
            else {
              inputs[key] = lastComponent.id
              inputValues[key] = lastComponent.data
            }
            likelihoods.push(scores[`${lastComponent.type} ${proc.spec}`] || 0)
          })
          const likelihood = mean(likelihoods)
          console.log(`${ctx.maxId} ${lastComponentProc.spec} => ${lastComponentProc.output.spec} => ${proc.spec}`)
          nextComponents.push({
            id: ++ctx.maxId,
            inputs,
            type: proc.spec,
            description: nodeDesc(proc, { inputs: inputValues }),
            likelihood,
          })
        })
    }
    if (nextComponents.length === 0) break
    else if (nextComponents.length === 1) {
      const [{ likelihood, ...nextComponent }] = nextComponents
      previousComponents.push(nextComponent)
      results.push({ role: 'assistant', content: '```component\n' + JSON.stringify(nextComponent) + '\n```', })
    }
    const totalLikelihood = nextComponents.reduce((totalLiklihood, component) => totalLiklihood + component.likelihood, 0.0)
    nextComponents.sort((a, b) => b.likelihood - a.likelihood)
    const likelyNextComponents = nextComponents.filter(nextComponent => (nextComponent.likelihood / totalLikelihood) >= 0.5)
    if (likelyNextComponents.length === 1) {
      const [{ likelihood: _, ...nextComponent }] = likelyNextComponents
      previousComponents.push(nextComponent)
      results.push({ role: 'assistant', content: '```component\n' + JSON.stringify(nextComponent) + '\n```', })
      continue
    } else {
      likelyNextComponents.forEach(nextComponent => {nextComponent.likelihood /= totalLikelihood})
      results.push({ role: 'assistant', content: '```suggestions\n' + JSON.stringify(nextComponents.slice(0, 3)) + '\n```', })
      break
    }
  }
  return results
}

export const ChatGPTEmbeddingCLF = API('/api/v1/chat/gpt-embedding-clf')
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
      newMessages.push(...(await findAnswers([...inputs.body.messages, ...newMessages])))
    } catch (e) {
      newMessages.push({ role: 'assistant', content: `Error: ${(e as any).toString()}` })
    }
    return newMessages
  })
  .build()
