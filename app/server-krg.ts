/**
 * The KRG for server-side code: the components (@/app/krg) along with the
 *  current user suggestions, read from the database on each call (suggestions can
 *  be added & removed at any time, by any instance).
 */
import db from '@/app/db'
import KRG from '@/core/KRG'
import { metanodes } from '@/components/server'
import { withSuggestions } from '@/app/extensions/suggestions'
import cache from '@/utils/global_cache'

const componentKRG = cache('server-component-krg', () => {
  const krg = new KRG()
  metanodes.forEach(metanode => krg.add(metanode))
  return krg
})

/** Complete component registry, without database-owned user suggestions. */
export function getComponentKRG() {
  return componentKRG
}

export async function getKRG() {
  const suggestions = await db.selectFrom('suggestion')
    .select(['name', 'inputs', 'output', 'description', 'user'])
    .execute()
  return withSuggestions(componentKRG, suggestions)
}
