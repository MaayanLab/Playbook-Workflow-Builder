import type KRG from '@/core/KRG'
import type { MetaNode } from '@/spec/metanode'

export function metanodeType(krg: KRG, metanode: MetaNode) {
  if (metanode.kind === 'data') {
    return krg.getNextProcess(metanode.spec).length === 0 ? 'Viewer' : 'Data'
  }
  const terminal = krg.getNextProcess(metanode.output.spec).length === 0
  if ('prompt' in metanode) return terminal ? 'Interactive Viewer' : 'Prompt'
  return terminal ? 'View Resolver' : 'Resolver'
}
