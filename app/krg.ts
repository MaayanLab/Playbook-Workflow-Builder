import KRG from '@/core/KRG'
import * as components from '@/components'
import { ensureArray } from '@/utils/array'
import { MetaNodeType } from '@/spec/metanode'

declare global {
  var krg: KRG | undefined
}

type ValuesOf<T> = T[keyof T]

const krg = global.krg || new KRG()
Object.values(components)
  .flatMap(component => ensureArray<ValuesOf<typeof components>>(component))
  .filter((component): component is MetaNodeType<unknown> => typeof component === 'object' && 'spec' in component)
  .forEach(metanode => krg.add(metanode))

global.krg = krg
export default krg
