import { decode_complete_process_input_refs, decode_complete_process_inputs, decode_complete_process_output, decode_complete_process_output_ref } from "./engine";
import extractCitations from "@/utils/citations";
import * as dict from '@/utils/dict'
import type { FPL } from "@/core/FPPRG"
import type KRG from "@/core/KRG"
import { Story } from "@/spec/metanode";

export type Metadata = {
  title?: string,
  description?: string,
}

export type Author = {
  name: string,
  affiliation?: string,
  email?: string,
  orcid?: string,
}

export async function fpl_expand(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null }) {
  const fullFPL = props.fpl.resolve()
  const processLookup = dict.init(
    await Promise.all(fullFPL.map(async (step, index) => {
      const metanode = props.krg.getProcessNode(step.process.type)
      let story: Story
      let inputs: Record<string, unknown> | undefined
      try { inputs = await decode_complete_process_inputs(props.krg, step.process) } catch (e) {}
      let input_refs: Record<string, string | string[]> | undefined
      try { input_refs = decode_complete_process_input_refs(props.krg, step) } catch (e) {}
      let output: unknown | undefined
      try { output = await decode_complete_process_output(props.krg, step.process) } catch (e) {}
      let output_ref = `\\figref{${step.id}}`
      story = metanode.story ? metanode.story({ inputs, output, input_refs, output_ref }) : {}
      return {
        key: step.process.id,
        value: {
          index,
          node: step.process,
          inputs, output,
          metanode,
          story,
        },
      }
    }))
  )
  const story = extractCitations(
    dict.items(processLookup).flatMap((proc) =>
      dict.items(proc.value.story).map((st) => ({
        text: st.value as string,
        tags: [proc.key as string, st.key as string],
      }))
    )
  )
  return {
    fullFPL,
    processLookup,
    story,
  }
}