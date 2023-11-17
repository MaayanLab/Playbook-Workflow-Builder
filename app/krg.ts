import KRG from '@/core/KRG'
import { metanodes } from '@/components'
import cache from '@/utils/global_cache'

const krg = cache('krg', () => new KRG())
metanodes.forEach(metanode => krg.add(metanode))

export default krg
