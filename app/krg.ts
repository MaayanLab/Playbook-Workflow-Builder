import KRG from '@/core/KRG'
import * as components from '@/components'
const krg = new KRG()
Object.values(components).forEach(krg.add)
export default krg