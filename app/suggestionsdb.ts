import { MemoryDatabase, PgDatabase, Database } from '@/app/extensions/suggestions'
declare global {
  var suggestionsdb: Database | undefined
}

let suggestionsdb: Database
if (process.env.DATABASE_URL) {
  if (!global.suggestionsdb) {
    global.suggestionsdb = suggestionsdb = new PgDatabase(process.env.DATABASE_URL)
  } else {
    suggestionsdb = global.suggestionsdb
  }
} else {
  if (!global.suggestionsdb) {
    global.suggestionsdb = suggestionsdb = new MemoryDatabase()
  } else {
    suggestionsdb = global.suggestionsdb
  }
}

export default suggestionsdb
