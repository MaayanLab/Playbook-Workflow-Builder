import handler from '@/utils/next-rest'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import { z } from 'zod'

const GPTRequest = z.object({
  model: z.string(),
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
    role: z.enum(['user', 'assistant']),
    content: z.string(),
  })),
})

const PROMPT = {
  model: 'gpt-3.5-turbo',
  temperature: 0.0,
  messages: [
    {
      role: 'system',
      content: `
Your task is to perform actions to interact with a user and ultimately answer their questions. You can perform the following actions:

FindComponent - Locate components relevant for the task at hand
ShowComponent - Show a component to the user as an answer to their question
Reply - Reply directly to the user to ask follow up questions

Use the following format:

Q: Question to be answered
T: Thought, you should always think about what to do
A: Action to take, should be one of ["FindComponent", "ShowComponent", "Reply"]
I: Input argument for the action
O: Observation, the result of the action
... (this Q/T/A/I/O will continue to repeat)

Begin!

Q: I'm looking for drugs which down regulate ACE2.
T: I should capture the gene the user provided.
A: Input
I: { "type": "Gene", "value": "ACE2" }
T: The user wants drugs, I should find a way to get drugs given a gene.
A: FindPaths
I: { "type": "Drugs" }
O: [{"name": "LINCSL1000ReverseSearchDashboard", "input": "Gene", "output": "Drugs", "description": "A dashboard for performing L1000 Reverse Search queries for a given gene"}]
`,
    },
  ],
  stop: [
    '\nO:', // don't try to make the observations yourself
    '\nQ:', // don't try to ask questions for the user
  ],
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
  const gptRes = await chatGPT({
    ...PROMPT,
    messages: [
      ...PROMPT.messages,
      ...body.messages.map(message =>
        message.role === 'user' ? {
          role: 'user', content: `Q: ${message.content}`,
        } : message
      ),
    ],
  })
  const gptResponse = {
    role: 'assistant',
    content: (
      gptRes.choices[0].message.content
        .replace(/^\n*/g, '')
        .replace(/\n*$/g, '')
    ),
  }
  res.status(200).send(JSON.stringify(gptResponse))
})
