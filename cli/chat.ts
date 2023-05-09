import { z } from 'zod'
import KRG from '@/core/KRG'
import * as components from '@/components'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { MetaNode, ProcessMetaNode } from '@/spec/metanode'
import readline from 'readline'

const krg = new KRG()
dict.values(components as any).flatMap(component => array.ensureArray(component))
  .filter((component): component is MetaNode => typeof component === 'object' && 'spec' in component)
  .forEach(metanode => krg.add(metanode))

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

async function openaiCompletions({ messages = [], model = 'gpt-3.5-turbo', temperature = 0.7, ...kwargs }: Partial<z.TypeOf<typeof OpenAICompletionsRequest>>) {
  console.log(`gpt> ${JSON.stringify(messages)}`)
    const req = await fetch(`https://api.openai.com/v1/chat/completions`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${process.env.OPENAI_API_KEY}`,
    },
    body: JSON.stringify({ messages, model, temperature, ...kwargs }),
  })
  return OpenAICompletionsResponse.parse(await req.json())
}

async function identifyInputSpec(question: string) {
  const res = await openaiCompletions({
    messages: [{
      role: 'user',
      content: [
        'Consider the following question:',
        question,
        'Consider choices below in the form [id]\t[description]:',
        ...krg.getNextProcess('').map(node =>
          `${node.spec}\t${node.output.meta.description}`
        ),
        '*Pluralality implies sets rather than terms.',
        'Respond with a JSON serialized mapping [id] to the applicable information from the question.',
      ].join('\n'),
    }]
  })
  return JSON.parse(strip(res.choices[0].message.content))
}

function nodeDesc(node: ProcessMetaNode, { input, output }: { input?: string, output?: string }) {
  if (node.story && Object.keys(node.inputs).length === 1) {
    try {
      node.story({
        inputs: { [Object.keys(node.inputs)[0] as any]: input || undefined },
        output,
      })
    } catch (e) {}
  }
  return node.meta.description
}

async function identifyNextStep(question: string, inputSpec: Record<string, string>) {
  const res = await openaiCompletions({
    messages: [{
      role: 'user',
      content: [
        'Consider the following question:',
        question,
        'Consider choices below in the form [id]\t[description]:',
        ...dict.items(inputSpec).flatMap(({ key, value }) =>
          krg.getNextProcess(krg.getProcessNode(key).output.spec).map(node =>
            `${node.spec}\t${nodeDesc(krg.getProcessNode(key), { output: value })}. ${nodeDesc(node, { input: value })}.`
          )
        ),
        'The most appropriate choice id to answer the question is:',
      ].join('\n'),
    }]
  })
  return strip(res.choices[0].message.content)
}

const std = readline.createInterface({ input: process.stdin, output: process.stdout })
std.question('User: ', async (question) => {
  const inputSpec = await identifyInputSpec(question)
  console.log({ inputSpec })
  if (inputSpec !== null) {
    const nextStep = await identifyNextStep(question, inputSpec)
    console.log({ nextStep })
  }
  std.close()
})
