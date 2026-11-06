import { createDatabase } from '@/db'
import cache from '@/utils/global_cache'

const backend = cache('db', () => createDatabase(process.env.DATABASE_URL))

/** The kysely query builder */
const db = backend.db
export default db

/** Event bus & job queue */
export const bus = backend.bus

export const dialect = backend.dialect

