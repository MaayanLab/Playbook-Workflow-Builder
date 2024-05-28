import type OpenAI from 'openai'
import { z } from 'zod'

const AgentC = z.object({
  message: z.string(),
  suggestions: z.array(z.object({ id: z.number(), value: z.string().optional() })),
})

const UserC = z.intersection(
  z.union([
    z.object({
      message: z.string(),
    }),
    z.object({
      step: z.object({
        id: z.number(),
        value: z.string().optional(),
      }),
    }),
  ]),
  z.object({
    choices: z.array(z.object({
      id: z.number(),
      name: z.string(),
      inputs: z.record(z.string(), z.object({ id: z.number() })),
      type: z.string(),
      label: z.string(),
      description: z.string(),
      story: z.string().optional(),
    })).optional(),
  }),
)

export function GPTAssistantMessageParse(messages: OpenAI.Beta.Threads.Messages.Message[]) {
  return messages.flatMap(msg =>
    msg.content.flatMap(content => {
      if (content.type === 'text') {
        if (msg.role === 'assistant') {
          const { data, error } = AgentC.safeParse(JSON.parse(content.text.value))
          if (data) return [{ role: 'assistant' as const, ...data }]
          else return [{ role: 'error' as const, error }]
        } else if (msg.role === 'user') {
          const { data, error } = UserC.safeParse(JSON.parse(content.text.value))
          if (data) return [{ role: 'user' as const, ...data }]
          else return [{ role: 'error' as const, error }]
        }
      }
      return [] as Array<
        ({ role: 'assistant' } & z.infer<typeof AgentC>)
        | ({ role: 'user' } & z.infer<typeof UserC>)
        | ({ role: 'error', error: any })
      >
    })
  )
}

export type AssistantParsedMessages = ReturnType<typeof GPTAssistantMessageParse>

export function AssembleState(messages: AssistantParsedMessages, { with_value = false } = {}) {
  let max_id = 0
  const all_nodes: Record<number, { id: number, name: string, value?: string, inputs: Record<string, { id: number }> }> = {}
  const workflow: { id: number, name: string, value?: string, inputs: Record<string, { id: number }> }[] = []
  messages
    .forEach(item => {
      if ('step' in item) {
        workflow.push(all_nodes[item.step.id])
        if (with_value && item.step.value) all_nodes[item.step.id].value = item.step.value
      }
      if ('choices' in item && item.choices) {
        item.choices.forEach(choice => {
          all_nodes[+choice.id] = choice
          max_id = Math.max(max_id, +choice.id)
        })
      }
    })
  return { all_nodes, workflow, max_id }
}
