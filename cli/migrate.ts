/**
 * Apply pending database migrations (postgres). sqlite databases are migrated automatically.
 *
 * Usage: DATABASE_URL=postgres://... npm run db:migrate
 */
import { createDatabase, migrateToLatest } from '@/db'

async function main() {
  const { db, dialect } = createDatabase(process.env.DATABASE_URL)
  if (dialect !== 'postgres') console.warn('DATABASE_URL is not postgres, sqlite databases are migrated automatically')
  try {
    await migrateToLatest(db)
    console.log('database is up to date')
  } finally {
    await db.destroy()
  }
  process.exit(0)
}

main().catch((error) => {
  console.error(error)
  process.exit(1)
})
