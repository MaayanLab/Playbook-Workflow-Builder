import KRG from '@/core/KRG'
import * as components from '@/components'
import { ensureArray } from '@/utils/array'
import { MetaNodeType } from '@/spec/metanode'
import cache from '@/utils/global_cache'

const krg = cache('krg', () => new KRG())
// @ts-ignore
Object.values(components).flatMap(component => ensureArray(component))
  .filter((component): component is MetaNodeType<unknown> => typeof component === 'object' && 'spec' in component)
  .forEach(metanode => krg.add(metanode))

export default krg
