/**
 * The KRG for server-side code: the components (@/app/krg) along with the
 *  current user suggestions, read from the database on each call (suggestions can
 *  be added & removed at any time, by any instance).
 */
import db from '@/app/db'
import krg from '@/app/krg'
import { withSuggestions } from '@/app/extensions/suggestions'

export async function getKRG() {
  const suggestions = await db.selectFrom('suggestion')
    .select(['name', 'inputs', 'output', 'description', 'user'])
    .execute()
  return withSuggestions(krg, suggestions)
}
