import handler from '@/utils/next-rest'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import { z } from 'zod'
import dedent from 'ts-dedent'

const GPTRequest = z.object({
  model: z.string(),
  temperature: z.number(),
  stop: z.array(z.string()),
  messages: z.array(z.object({
    role: z.string(),
    content: z.string(),
  })),
})

const GPTResponse = z.object({
  choices: z.array(z.object({
    message: z.object({
      role: z.string(),
      content: z.string(),
    })
  }))
})

async function chatGPT(body: z.TypeOf<typeof GPTRequest>) {
  const gptReq = await fetch(`https://api.openai.com/v1/chat/completions`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${process.env.OPENAI_API_KEY}`,
    },
    body: JSON.stringify(body),
  })
  const gptRes = GPTResponse.parse(await gptReq.json())
  return gptRes
}

/**
 * 
 * Q: I'm looking for drugs which down regulate ACE2.
 * T: I need to find a component to help the user answer their question.
 * A: FindComponents
 * I: { "input": [{ "type": "Gene", "data": "ACE2" }], "output": [{ "type": "Drug" }] }
 * O: [{"name": "LINCSL1000ReverseSearchDashboard", "description": "A dashboard for performing L1000 Reverse Search queries for a given gene"}]
 * T: I can reply to the user.
 * R: LINCSL1000ReverseSearchDashboard
 */

const BodyType = z.object({
  messages: z.array(z.object({
    role: z.enum(['system', 'user', 'assistant']),
    content: z.string(),
  })),
})

function find(content: string) {
  return [
    { id: 'Gene', description: 'A gene' },
    { id: 'Drug', description: 'A drug' },
    { id: 'Metabolite', description: 'A metabolite' },
    { id: 'Variant', description: 'A genomic variant' },
  ]
}

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
  if (lastMessage.role !== 'user') throw new Error('Unexpected input')
  const content = [
    dedent`
      You are a workflow builder, you're helping a user construct a workflow.
      Your first goal is to figure out what the user is looking for and what they can provide.

      Use the following format:
      Q: the user's message
      O: possible types and their descriptions
      A: confirm with the user, of the possible types, what they're looking for and what they can provide
      ... Repeat Q/O/A until you're certain of the types.
      F: Looking For: {all,type,ids} Can Provide: {all,type,ids}

      Begin!
    `,
    ...body.messages.map(message => {
      if (message.role === 'user') return `Q: ${message.content}`
      else if (message.role === 'assistant') return `A: ${message.content}`
      else return message.content
    }),
    `O: ${JSON.stringify(find(lastMessage.content))}`
  ].join('\n')
  const gptRes = await chatGPT({
    model: 'gpt-3.5-turbo',
    temperature: 0.0,
    messages: [
      {
        role: 'user',
        content,
      },
    ],
    stop: [
      '\nQ:',
      '\nO:',
    ],
  })
  const gptResponse = {
    nodes: {},
    messages: [{
      role: 'assistant',
      content: gptRes.choices[0].message.content
        .replace(/^\n*/g, '')
        .replace(/\n*$/g, '')
    }]
  }
  res.status(200).json(gptResponse)
})
