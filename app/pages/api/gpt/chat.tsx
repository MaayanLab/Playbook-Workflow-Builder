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

async function identifyInputs(question: string) {
  const res = await openaiCompletions({
    messages: [{
      role: 'user',
      content: [
        'Consider the following user query:',
        question,
        'Consider the choices below:',
        JSON.stringify(krg.getNextProcess('').map((node, i) => ({ id: node.spec, description: node.output.meta.description }))),
        '*Pluralality implies sets rather than terms.',
        'Respond with a precise JSON serialized mapping choice id to the applicable information from the query.',
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
        inputs: { [Object.keys(node.inputs)[0] as any]: input || `[${Object.keys(node.inputs)[0]}]` },
        output,
      })
    } catch (e) {}
  }
  return node.meta.description
}

async function identifyWorkflow(question: string, inputs: Record<string, string>) {
  const workflows = dict.items(inputs).flatMap(({ key, value }) => {
    const inputNode = krg.getProcessNode(key)
    if (!inputNode) return []
    return krg.getNextProcess(inputNode.output.spec).map(node => ({
      workflow: [
        { spec: key, data: value },
        { spec: node.spec },
      ],
      description: `${nodeDesc(inputNode, { output: value })}. ${nodeDesc(node, { input: value })}.`,
    }))
  })
  const res = await openaiCompletions({
    messages: [{
      role: 'user',
      content: [
        'Consider the following user query:',
        question,
        'Consider the workflows below:',
        JSON.stringify(workflows.map((choice, i) => ({ id: i, ...choice }))),
        'Response should be a JSON serialization which satisfies the typescript type:',
        `{ /* most appropriate workflow id satisfying the query */ "workflow": number } | { /* ask the user if no choices are applicable */ "question": string }`,
      ].join('\n'),
    }]
  })
  const content = strip(res.choices[0].message.content)
  try {
    const m = /\{.+?\}/gms.exec(content)
    const contentParsed = JSON.parse(m ? m[0] : '')
    if (typeof contentParsed.workflow === 'number') {
      const { workflow = undefined } = workflows[contentParsed.workflow] || {}
      if (workflow) return { success: true, content: '```workflow\n' + JSON.stringify({ step: 0, workflow }) + '\n```' }
    }
    return { success: false, content: contentParsed.question }
  } catch (e) {}
  return { success: false, content }
}

async function identifyFollowupWorkflow(messages: { role: "user" | "system" | "assistant", content: string }[]) {
  console.debug(messages)
  const [lastWorkflowMessage] = messages.filter(({ role, content }) => role === 'assistant' && content.includes('```workflow\n')).slice(-1)
  const m = /```workflow\n(.+)\n```/gms.exec(lastWorkflowMessage.content)
  const lastWorkflow = m ? JSON.parse(m[1]) : null
  const lastProcessNode = krg.getProcessNode(lastWorkflow.workflow[lastWorkflow.workflow.length-1].spec)
  const workflows = krg.getNextProcess(lastProcessNode.output.spec).map(node => {
    const workflow = [...lastWorkflow.workflow, { spec: node.spec }]
    return {
      workflow,
      description: workflow.map(
        ({ spec }, i) =>
          nodeDesc(krg.getProcessNode(spec), {
            input: workflow[i].data,
            output: i > 0 ? workflow[i-1].data : undefined
          })).join('. ') + '.',
    }
  })
  if (workflows.length === 0) return { success: false, content: '' }
  const res = await openaiCompletions({
    messages: [
      ...messages.slice(0, -1),
      {
        role: 'user',
        content: [
          'Also considering our previous conversation, consider the following query:',
          messages[messages.length-1],
          'Consider the following potential next workflows:',
          JSON.stringify(workflows.map((workflow, i) => ({ id: i, ...workflow }))),
          'Response should be a JSON serialization which satisfies the typescript type:',
          `{ /* most appropriate workflow id satisfying the query */ "workflow": number } | { /* ask the user if no choices are applicable */ "question": string }`,
        ].join('\n'),
      },
    ]
  })
  const content = strip(res.choices[0].message.content)
  try {
    const m = /\{.+?\}/gms.exec(content)
    const contentParsed = JSON.parse(m ? m[0] : '')
    if (typeof contentParsed.workflow === 'number') {
      const { workflow = undefined } = workflows[contentParsed.workflow] || {}
      if (workflow) return { success: true, content: '```workflow\n' + JSON.stringify({ step: lastWorkflow.workflow.length, workflow }) + '\n```' }
    }
    return { success: false, content: contentParsed.question }
  } catch (e) {}
  return { success: false, content }
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
  if (body.messages.length === 1) {
    const inputs = await identifyInputs(lastMessage.content)
    let step = inputs.success ? await identifyWorkflow(lastMessage.content, inputs.content) : inputs
    let previousStep: typeof step | undefined
    // while (step.success) {
    //   previousStep = step
    //   body.messages.splice(body.messages.length-1, 0, { role: 'assistant', content: step.content })
    //   step = await identifyFollowupWorkflow(body.messages)
    // }
    res.status(200).json({ role: 'assistant', content: [previousStep?.content, step.content].filter(c => c!==undefined).join('\n') })
  } else {
    let step = await identifyFollowupWorkflow(body.messages)
    let previousStep: typeof step | undefined
    // while (step.success) {
    //   previousStep = step
    //   body.messages.splice(body.messages.length-1, 0, { role: 'assistant', content: step.content })
    //   step = await identifyFollowupWorkflow(body.messages)
    // }
    res.status(200).json({ role: 'assistant', content: [previousStep?.content, step.content].filter(c => c!==undefined).join('\n') })
  }
})
