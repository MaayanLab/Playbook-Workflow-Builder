import type KRG from '@/core/KRG'
import type { Metapath } from '@/app/fragments/graph/types'
import { useMetapathInputs, useMetapathOutput } from './metapath'

export default function Story({ krg, head }: { krg: KRG, head: Metapath }) {
  const processNode = krg.getProcessNode(head.process.type)
  const { data: inputs } = useMetapathInputs(krg, head)
  const { data: { output } } = useMetapathOutput(krg, head)
  if (!inputs || !processNode.story) return null
  try {
    return (
      <>
        {processNode.story({
          inputs,
          output,
        })}&nbsp;
      </>
    )
  } catch (e) {
    return null
  }
}
