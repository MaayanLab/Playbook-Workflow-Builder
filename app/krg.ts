import KRG from '@/core/KRG'
import * as components from '@/components'
import { ensureArray } from '@/utils/array'
import { MetaNodeType } from '@/spec/metanode'
import cache from '@/utils/global_cache'

type ValuesOf<T> = T[keyof T]

const krg = cache('krg', () => new KRG())

Object.values(components)
  .flatMap(component => ensureArray<ValuesOf<typeof components>>(component))
  .filter((component): component is MetaNodeType<unknown> => typeof component === 'object' && 'spec' in component)
  .forEach(metanode => krg.add(metanode))

export default krg
