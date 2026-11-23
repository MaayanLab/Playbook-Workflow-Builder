/**
 * How user_playbooks are indexed for search (see db/queries.ts, db/migrations/0010_user_playbook_search.ts).
 *  Free of db imports so that migrations can use it.
 */
import { tsvector } from '../utils/tsvector'

export type UserPlaybookSearchable = { title?: string | null, description?: string | null, inputs?: string | null, outputs?: string | null }

/** Data types are stored joined, e.g. `Term[Gene], Set[Gene]` */
export function splitTypes(types: string | null | undefined) {
  return types ? types.split(', ').filter(Boolean) : []
}

/** The rows indexing a user_playbook for search: the trigrams of its text & its input/output types */
export function userPlaybookIndexRows(id: string, record: UserPlaybookSearchable) {
  const text = [record.title || '', record.description || '', record.inputs || '', record.outputs || ''].join(' ')
  return {
    trigrams: [...tsvector(text)].map(trigram => ({ user_playbook: id, trigram })),
    types: [
      ...new Set(splitTypes(record.inputs).map(spec => `input\0${spec}`)),
      ...new Set(splitTypes(record.outputs).map(spec => `output\0${spec}`)),
    ].map(key => {
      const [kind, spec] = key.split('\0') as ['input' | 'output', string]
      return { user_playbook: id, kind, spec }
    }),
  }
}
