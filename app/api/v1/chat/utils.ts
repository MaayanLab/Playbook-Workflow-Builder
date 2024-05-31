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

const zJson = z.string().transform((content, ctx) => {
  try {
    const m = /({.+})/sg.exec(content)
    if (!m) return content
    return JSON.parse(m[1])
  } catch (e) {
    ctx.addIssue({
      code: z.ZodIssueCode.invalid_string,
      validation: 'json' as any,
      message: (e as Error).message,
    })
    return z.NEVER
  }
})

function zodCodecTransform<Output = any, Def extends z.ZodTypeDef = z.ZodTypeDef, Input = Output>(zodCodec: z.ZodType<Output, Def, Input>) {
  return (content: Input, ctx: z.RefinementCtx): Output => {
    const { data, error } = zodCodec.safeParse(content)
    if (error) {
      error.issues.forEach(issue => ctx.addIssue(issue))
      return z.NEVER
    }
    return data
  }
}

export function GPTAssistantMessageParse(messages: { id: string, role: string, content: string }[]) {
  return messages.flatMap((msg) => {
    if (msg.role === 'assistant') {
      const { data, error } = zJson.transform(zodCodecTransform(AgentC)).safeParse(msg.content)
      if (data) return [{ id: msg.id, role: 'assistant' as const, ...data }]
      else return [{ id: msg.id, role: 'error' as const, message: error.toString() }]
    } else if (msg.role === 'user') {
      const { data, error } = zJson.transform(zodCodecTransform(UserC)).safeParse(msg.content)
      if (data) return [{ id: msg.id, role: 'user' as const, ...data }]
      else return [{ id: msg.id, role: 'error' as const, message: error.toString() }]
    }
    return [] as Array<
      ({ id: string, role: 'assistant' } & z.infer<typeof AgentC>)
      | ({ id: string, role: 'user' } & z.infer<typeof UserC>)
      | ({ id: string, role: 'error', message: string })
    >
  })
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
