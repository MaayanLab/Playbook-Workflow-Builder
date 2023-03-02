import type KRG from '@/core/KRG'
import type { FPL } from '@/core/FPPRG'
import * as dict from '@/utils/dict'
import pluralize from 'pluralize'

type Metapath = ReturnType<FPL['toJSON']>

function augment(s: string, capitalize = false) {
  const m = /^([A-Z]?)(.+)(\.?)$/.exec(s)
  if (m) return `${capitalize ? m[1].toUpperCase() : m[1].toLowerCase()}${m[2]}`
  else return s
}

function sentence_comma_separate(s: string[]) {
  if (s.length === 1) return `the ${s}`
  else if (s.length === 2) return s.join(' and ')
  else if (s.length >= 3) return `${s.slice(0, -1).join(', ')}, and ${s[s.length-1]}`
}

export default function FPL2Text(krg: KRG, metapath: Array<Metapath>) {
  return metapath.map((head, index) => {
    const proc = krg.getProcessNode(head.process.type)
    if (dict.isEmpty(proc.inputs)) {
      return `${augment(proc.meta.description, true)}.`
    } else {
      return `Using ${sentence_comma_separate(
        dict.values(proc.inputs).map(input => {
          if (Array.isArray(input)) {
            return `several ${pluralize(augment(input[0].meta.description))}`
          } else {
            return augment(input.meta.description)
          }
        })
      )}, ${augment(proc.meta.description)} to produce ${augment(proc.output.meta.description)}.`
    }
  }).join(' ')
}
