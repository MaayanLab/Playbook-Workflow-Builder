import KRG from '@/core/KRG'
import * as components from '@/components'
import { ensureArray } from '@/utils/array'
import { MetaNodeType } from '@/spec/metanode'

declare global {
  var krg: KRG | undefined
}

const krg = global.krg || new KRG()
// @ts-ignore
Object.values(components).flatMap(component => ensureArray(component))
  .filter((component): component is MetaNodeType<unknown> => typeof component === 'object' && 'spec' in component)
  .forEach(metanode => krg.add(metanode))

global.krg = krg
export default krg
