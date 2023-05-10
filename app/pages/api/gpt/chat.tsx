import krg from '@/app/krg'
import handler from '@/utils/next-rest'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
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

async function identifyInputSpec(question: string) {
  const res = await openaiCompletions({
    messages: [{
      role: 'user',
      content: [
        'Consider the following question:',
        question,
        'Consider the choices below:',
        JSON.stringify(krg.getNextProcess('').map((node, i) => ({ id: node.spec, description: node.output.meta.description }))),
        '*Pluralality implies sets rather than terms.',
        'Respond with a precise JSON serialized mapping choice id to the applicable information from the question.',
      ].join('\n'),
    }]
  })
  try {
    return { success: true, content: JSON.parse(strip(res.choices[0].message.content)) }
  } catch (e) {
    return { success: false, content: strip(res.choices[0].message.content) }
  }
}

function nodeDesc(node: ProcessMetaNode, { input, output }: { input?: string, output?: string }) {
  if (node.story && Object.keys(node.inputs).length === 1) {
    try {
      node.story({
        inputs: { [Object.keys(node.inputs)[0] as any]: input || '' },
        output,
      })
    } catch (e) {}
  }
  return node.meta.description
}

async function identifyWorkflow(question: string, inputSpec: Record<string, string>) {
  const choices = dict.items(inputSpec).flatMap(({ key, value }) => {
    const inputNode = krg.getProcessNode(key)
    if (!inputNode) return []
    return krg.getNextProcess(inputNode.output.spec).map(node => ({
      input: key, inputData: value, output: node.spec,
      description: `${nodeDesc(inputNode, { output: value })}. ${nodeDesc(node, { input: value })}.`,
    }))
  })
  const res = await openaiCompletions({
    messages: [{
      role: 'user',
      content: [
        'Consider the following question:',
        question,
        'Consider the choices below:',
        JSON.stringify(choices.map((choice, i) => ({ id: i, description: choice.description }))),
        'Respond only with the most appropriate choice id for addressing the question.',
      ].join('\n'),
    }]
  })
  const choice = choices[Number(strip(res.choices[0].message.content))]
  if (!choice) return { success: false, content: strip(res.choices[0].message.content) }
  else return { success: true, content: '```workflow\n' + JSON.stringify([{ spec: choice.input, data: choice.inputData }, { spec: choice.output }]) + '\n```' }
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
  const inputs = await identifyInputSpec(lastMessage.content)
  const { content } = inputs.success ? await identifyWorkflow(lastMessage.content, inputs.content) : inputs
  res.status(200).json({ role: 'assistant', content })
})
