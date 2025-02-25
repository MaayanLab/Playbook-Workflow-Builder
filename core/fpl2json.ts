import type { FPL } from "./FPPRG";
import * as dict from '@/utils/dict'

export default async function fpl2json(props: { fpl: FPL }) {
  const data: Record<string, any> = {}
  const workflow = []
  for (const el of props.fpl.resolve()) {
    const output = await el.process.output()
    workflow.push({
      id: el.process.id,
      type: el.process.type,
      inputs: !dict.isEmpty(el.process.inputs) ? dict.init(dict.items(el.process.inputs).map(({ key, value }) => ({ key, value: { id: value.id } }))) : undefined,
      data: el.process.data ? { id: el.process.data.id } : undefined,
      output: output ? { id: output.id } : undefined
    })
    if (el.process.data) {
      const { id: _, ...data_node } = el.process.data.toJSON()
      data[el.process.data.id] = data_node
    }
    if (output) {
      const { id: _, ...data_node } = output.toJSON()
      data[output.id] = data_node
    }
  }
  return {
    data,
    workflow,
    metadata: props.fpl.playbook_metadata ? props.fpl.playbook_metadata.toJSON() : undefined
  }
}
