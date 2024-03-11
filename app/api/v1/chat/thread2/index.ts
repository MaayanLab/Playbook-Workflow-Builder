import { API } from '@/spec/api'
import { z } from 'zod'
import krg from '@/app/krg'
import OpenAI from 'openai'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { ResponseCodedError, UnauthorizedError } from '@/spec/error'
import cache from '@/utils/global_cache'
import dedent from 'ts-dedent'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { ProcessMetaNode } from '@/spec/metanode'

const openai = cache('openai', () => new OpenAI({ apiKey: process.env.OPENAI_API_KEY }))
const assistant = cache('openai-assistant2', async () => {
  if (!process.env.OPENAI_ASSISTANT_ID) {
    const assistant = await openai.beta.assistants.create({
      name: `playbook-gpt`,
      instructions: dedent`
        You are a bot that helps the user construct workflows. You do not answer questions directly, instead you help choose from the applicable workflow steps in an effort to construct a workflow that answer the user's question.
        Queries will be JSON serialized, and of the form (note that steps in the current workflow are linked to the next step by id):
        User:
        {
          "message": "I'd like to get the expression of ACE2.",
          "choices": [{ "id": 1, "type": "prompt", "name": "Input[Gene]", ...]
        }
        Assistant:
        {
          "message": "In that case, let's start with the gene, ACE2.",
          "suggestions": [{ "id": 1, "value": "ACE2" }, ...]
        }
        User:
        {
          "step": { "id": 1 },
          "choices": [{ "id": 2, "type": "resolver", "name": "GetTissueExpressionOfGene", "inputs": { "gene": 1 }, ... }, ...]
        }
        Assistant:
        {
          "message": "You can get tissue expression of ACE2 in healthy human tissue from GTEx or from GEO using ARCHS4.",
          "suggestions": [{ "id": 2 }, ...]
        }
        Your response should be JSON serialized, choosing from these or choices in past messages, suggestions should be ordered by most likely.
        When the choice type is prompt, value must be specified based on the user search term.
      `,
      model: "gpt-4-turbo-preview",
    })
    console.log(assistant.id)
    return assistant
  } else {
    return await openai.beta.assistants.retrieve(process.env.OPENAI_ASSISTANT_ID)
  }
})

function TryJSONParse(obj: string) {
  try {
    return JSON.parse(obj)
  } catch (e) {
    return obj
  }
}

function AssembleState(messages: (string | undefined)[]) {
  let max_id = 0
  const all_nodes: Record<number, { id: number, name: string, inputs: Record<string, { id: number }> }> = {}
  const workflow: { id: number, name: string, inputs: Record<string, { id: number }> }[] = []
  messages
    .filter((msg): msg is string => !!msg)
    .map(TryJSONParse)
    .map(msg => z.union([
      z.string().transform(el => ({ message: el })),
      z.object({ message: z.string(), choices: z.array(z.object({ id: z.number(), name: z.string(), inputs: z.record(z.string(), z.object({ id: z.number() })) })) }),
      z.object({ step: z.object({ id: z.number() }), choices: z.array(z.object({ id: z.number(), name: z.string(), inputs: z.record(z.string(), z.object({ id: z.number() })) })) }),
    ]).parse(msg))
    .forEach(item => {
      if ('step' in item) {
        workflow.push(all_nodes[+item.step.id])
      }
      if ('choices' in item) {
        item.choices.forEach(choice => {
          all_nodes[+choice.id] = choice
          max_id = Math.max(max_id, +choice.id)
        })
      }
    })
  return { all_nodes, workflow, max_id }
}

function count(start: number) {
  const ctx = { id: start }
  return {
    next: () => ctx.id++
  }
}

export const GPTAssistantThread2Create = API.post('/api/v1/chat/threads2')
  .query(z.object({}))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const thread = await openai.beta.threads.create()
    return thread.id
  })
  .build()

export const GPTAssistantThread2Message = API.post('/api/v1/chat/threads2/[thread_id]/messages')
  .query(z.object({ thread_id: z.string() }))
  .body(z.union([
    z.object({
      message: z.string(),
    }),
    z.object({
      step: z.object({ id: z.number() }),
    }),
  ]))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const currentMessages = await openai.beta.threads.messages.list(inputs.query.thread_id, { order: 'asc' })
    const { all_nodes, workflow, max_id } = AssembleState(
      currentMessages.data
        .filter(msg => msg.role === 'user')
        .flatMap(msg => msg.content)
        .map(msgContent => msgContent.type === 'text' ? msgContent.text.value : undefined)
    )
    if ('step' in inputs.body) workflow.push(all_nodes[inputs.body.step.id])
    const last = workflow.length > 0 ? krg.getProcessNode(workflow[workflow.length-1].name).output.spec : ''
    const id = count(max_id+1)
    const selections: Record<string, { process: typeof workflow[0], processNode: ProcessMetaNode }> = {}
    workflow.forEach(item => {
      if (item === undefined) return
      // add this to the selections
      selections[item.id] = { process: item, processNode: krg.getProcessNode(item.name) }
      // if a selection previously registered is a parent of this selection, remove it from selections
      dict.values(item.inputs).forEach(k => {
        if (k.id in selections) delete selections[k.id]
      })
    })
    const choices = krg.getNextProcess(last)
      .filter(item => item.meta.hidden !== true)
      .filter(item => array.all(
        dict.values(item.inputs).map((value) => {
          if (Array.isArray(value)) {
            return dict.values(selections).filter(selection => selection.processNode.output.spec === value[0].spec).length > 1
          } else {
            return dict.values(selections).filter(selection => selection.processNode.output.spec === value.spec).length >= 1
          }
        })
      ))
      .map(item => {
        const inputs: Record<string, { id: number }> = {}
        dict.items(item.inputs).forEach(({ key: arg, value: input }) => {
          if (Array.isArray(input)) {
            dict.values(selections)
              .filter(selection => selection.processNode.output.spec === input[0].spec)
              .forEach((selection, i) => {
                inputs[`${arg}:${i}`] = { id: selection.process.id }
              })
          } else {
            const head = { process: workflow[workflow.length-1] }
            const relevantSelections = dict.filter(selections, ({ value: selection }) => selection.processNode.output.spec === input.spec)
            const selection = head.process.id in relevantSelections ? head : array.ensureOne(dict.values(relevantSelections))
            inputs[arg] = { id: selection.process.id }
          }
        })
        return {
          id: id.next(), name: item.spec, inputs,
          type: 'prompt' in item ? 'prompt' : 'resolver', label: item.meta.label,
          description: item.meta.description, story: item.story({}),
        }
      })
    const userMessage = await openai.beta.threads.messages.create(inputs.query.thread_id, {
      role: 'user',
      content: JSON.stringify({
        ...inputs.body,
        choices,
      }),
    })
    console.log(JSON.stringify({ userMessage }))

    let run = await openai.beta.threads.runs.create(inputs.query.thread_id, { assistant_id: (await assistant).id })
    while (true) {
      await new Promise<void>((resolve, reject) => {setTimeout(() => {resolve()}, 500)})
      run = await openai.beta.threads.runs.retrieve(inputs.query.thread_id, run.id)
      if (run.status === 'completed') {
        // send all new messages since the user's message to the user
        const newMessages = await openai.beta.threads.messages.list(inputs.query.thread_id, { after: userMessage.id, order: 'asc' })
        return { data: newMessages.data }
      } else if (run.status === 'in_progress' || run.status === 'queued') {
        continue
      } else {
        throw new ResponseCodedError(500, 'An OpenAI error occurred, try again later')
      }
    }
  })
  .build()

export const GPTAssistantThread2MessagesList = API.get('/api/v1/chat/threads2/[thread_id]/messages')
  .query(z.object({ thread_id: z.string() }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    const messages = await openai.beta.threads.messages.list(inputs.query.thread_id, { order: 'asc' })
    return messages.data.map(datum => ({
      role: datum.role,
      content: datum.content,
    }))
  })
  .build()

export const GPTAssistantThread2Delete = API.post('/api/v1/chat/threads2/[thread_id]/delete')
  .query(z.object({ thread_id: z.string() }))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || !process.env.OPENAI_API_KEY) throw new UnauthorizedError()
    if (req.method === 'HEAD') {
      res.status(200).end()
      return
    }
    await openai.beta.threads.del(inputs.query.thread_id)
  })
  .build()
