/**
 * Migrations are registered explicitly (rather than read from disk)
 *  so that they work when bundled by next.js.
 *
 * Add new migrations here with a lexicographically increasing key.
 */
import type { Migration, MigrationProvider } from 'kysely/migration'
import * as m0000_baseline from './0000_baseline'
import * as m0001_thread_openai_thread_nullable from './0001_thread_openai_thread_nullable'
import * as m0002_system_user from './0002_system_user'
import * as m0003_pgboss_v12 from './0003_pgboss_v12'
import * as m0004_proxy_session_user from './0004_proxy_session_user'
import * as m0005_thread_message_response_item from './0005_thread_message_response_item'
import * as m0006_better_auth from './0006_better_auth'
import * as m0007_fpl_export from './0007_fpl_export'
import * as m0008_api_keys from './0008_api_keys'

const migrations: Record<string, Migration> = {
  '0000_baseline': m0000_baseline,
  '0001_thread_openai_thread_nullable': m0001_thread_openai_thread_nullable,
  '0002_system_user': m0002_system_user,
  '0003_pgboss_v12': m0003_pgboss_v12,
  '0004_proxy_session_user': m0004_proxy_session_user,
  '0005_thread_message_response_item': m0005_thread_message_response_item,
  '0006_better_auth': m0006_better_auth,
  '0007_fpl_export': m0007_fpl_export,
  '0008_api_keys': m0008_api_keys,
}

export const migrationProvider: MigrationProvider = {
  getMigrations: async () => migrations,
}
