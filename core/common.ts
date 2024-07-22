import { decode_complete_process_inputs, decode_complete_process_output } from "./engine";
import extractCitations from "@/utils/citations";
import * as dict from '@/utils/dict'
import type { FPL } from "@/core/FPPRG"
import type KRG from "@/core/KRG"
import { Story } from "@/spec/metanode";
import Citable from "@/utils/citations";

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
      let story: Citable | undefined
      let inputs: Record<string, unknown> | undefined
      let output: unknown | undefined
      try { inputs = await decode_complete_process_inputs(props.krg, step.process) } catch (e) {}
      try { output = await decode_complete_process_output(props.krg, step.process) } catch (e) {}
      if (metanode.story) {
        story = Citable.concat(dict.items(metanode.story({ inputs, output })).map(({ key, value}) => Citable.ensure(value).section(key)))
      }
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
  const story = Citable.concat(dict.values(processLookup).map(({ story }) => story))
  return {
    fullFPL,
    processLookup,
    story,
  }
}