import db from '@/app/db'

export async function workflowExists(id?: string) {
  if (!id || id === 'start') return true
  return await db.selectFrom('fpl').select('id').where('id', '=', id).executeTakeFirst() !== undefined
}

export async function workflowNodeExists(graphId: string, nodeId?: string) {
  if (!nodeId || nodeId === 'start') return await workflowExists(graphId)
  const ids = new Set(['start'])
  let id: string | null = graphId
  while (id) {
    const row: { id: string, parent: string | null } | undefined = await db
      .selectFrom('fpl')
      .select(['id', 'parent'])
      .where('id', '=', id)
      .executeTakeFirst()
    if (!row) return false
    ids.add(row.id)
    id = row.parent
  }
  return nodeId.split(',').every(id => ids.has(id))
}
