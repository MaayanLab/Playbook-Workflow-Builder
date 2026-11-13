import db from '@/app/db'
import fpprg from '@/app/fpprg'

export async function workflowExists(id?: string) {
  if (!id || id === 'start') return true
  return await db.selectFrom('fpl').select('id').where('id', '=', id).executeTakeFirst() !== undefined
}

/**
 * Whether `nodeId` (comma separated) are steps of workflow `graphId`.
 *
 * The steps are those of the workflow as it's served (fpprg.getFPL): a workflow stored with ids
 *  that its content no longer hashes to (e.g. one registered by an older version, like the demo
 *  playbooks) is re-addressed when it's read, so its steps have new ids, which the pages link to.
 *  The ids of the steps as stored are accepted too, for links made with them.
 */
export async function workflowNodeExists(graphId: string, nodeId?: string) {
  if (!nodeId || nodeId === 'start') return await workflowExists(graphId)
  const fpl = await fpprg.getFPL(graphId)
  if (!fpl) return false
  const ids = new Set(['start', ...fpl.resolve().map(step => step.id)])
  let id: string | null = graphId
  while (id) {
    const row: { id: string, parent: string | null } | undefined = await db
      .selectFrom('fpl')
      .select(['id', 'parent'])
      .where('id', '=', id)
      .executeTakeFirst()
    if (!row) break
    ids.add(row.id)
    id = row.parent
  }
  return nodeId.split(',').every(id => ids.has(id))
}
