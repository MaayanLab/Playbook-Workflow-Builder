import handler from '@/utils/next-rest'
import { getServerSessionWithId } from "@/app/extensions/next-auth/helpers"
import { UnauthorizedError, UnsupportedMethodError } from '@/spec/error'
import { z } from 'zod'

const BodyType = z.object({
  description: z.string(),
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
  const gptReq = await fetch(`https://api.openai.com/v1/chat/completions`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${process.env.OPENAI_API_KEY}`,
    },
    body: JSON.stringify({
      "model": "gpt-3.5-turbo",
      "messages": [{
        "role": "user",
        "content": `Restructure this text so it is easier to read, do not change or omit any information:\n${body.description}`
      }],
      "temperature": 0.7
    }),
  })
  const gptRes = z.object({
    choices: z.array(z.object({
      message: z.object({
        role: z.string(),
        content: z.string(),
      })
    })),
  }).parse(await gptReq.json())
  const gptResponse = (
    `${gptRes.choices[0].message.content
      .replace(/^\n*/g, '')
      .replace(/\n*$/g, '')} (ChatGPT)`
  )
  res.status(200).send(JSON.stringify(gptResponse))
})
