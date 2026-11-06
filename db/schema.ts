/**
 * Kysely table definitions shared by the postgres (production) and
 *  sqlite (development / ephemeral) backends.
 *
 * Conventions which keep both dialects returning identical values:
 *  - ids are uuid strings, always generated application-side (see `uuid()` in @/db)
 *  - timestamps are ISO-8601 strings (the pg pool parses `timestamp` into ISO strings)
 *  - bigint columns are returned as javascript numbers (the pg pool parses int8)
 *  - booleans are returned as booleans (the sqlite driver wrapper converts 0/1)
 */
import type { ColumnType, Insertable, Selectable, Updateable } from 'kysely'

/** A timestamp column with a database-side default of now() */
export type CreatedTimestamp = ColumnType<string, string | undefined, string>
/** A counter column with a database-side default of 0 */
export type Counter = ColumnType<number, number | undefined, number>
/** A nullable column which may be omitted on insert */
export type Nullable<T> = ColumnType<T | null, T | null | undefined, T | null>

// accounts

export interface UserTable {
  id: string
  name: Nullable<string>
  email: string
  emailVerified: Nullable<string>
  image: Nullable<string>
  affiliation: Nullable<string>
  created: CreatedTimestamp
}

export interface AccountTable {
  id: string
  userId: string
  type: 'oauth' | 'email' | 'credentials' | (string & {})
  provider: string
  providerAccountId: string
  refresh_token: Nullable<string>
  access_token: Nullable<string>
  expires_at: Nullable<number>
  token_type: Nullable<string>
  scope: Nullable<string>
  id_token: Nullable<string>
  session_state: Nullable<string>
  oauth_token_secret: Nullable<string>
  oauth_token: Nullable<string>
}

export interface SessionTable {
  id: string
  expires: Nullable<string>
  sessionToken: string
  userId: string
}

export interface VerificationTokenTable {
  id: string
  identifier: Nullable<string>
  token: Nullable<string>
  expires: Nullable<string>
}

export interface UserIntegrationsTable {
  id: string
  cavatica_api_key: Nullable<string>
  cavatica_default_project: Nullable<string>
}

export interface ProxySessionTable {
  id: string
  /** the user who launched the session */
  user: Nullable<string>
  state: Nullable<string>
  run_id: Nullable<string>
  created: CreatedTimestamp
}

export interface SuggestionTable {
  id: string
  user: string
  name: string
  inputs: string
  output: string
  description: string
  created: CreatedTimestamp
}

// fpprg

export interface DataTable {
  id: string
  type: string
  value: string
  created: CreatedTimestamp
  hits: Counter
}

export interface ProcessTable {
  id: string
  type: string
  data: Nullable<string>
  timestamp: Nullable<string>
  created: CreatedTimestamp
  hits: Counter
}

export interface ProcessInputTable {
  id: string
  key: string
  value: string
}

export interface ResolvedTable {
  id: string
  data: Nullable<string>
  created: CreatedTimestamp
}

export interface CellMetadataTable {
  id: string
  label: Nullable<string>
  description: Nullable<string>
  process_visible: Nullable<boolean>
  data_visible: Nullable<boolean>
  created: CreatedTimestamp
}

export interface PlaybookMetadataTable {
  id: string
  title: Nullable<string>
  description: Nullable<string>
  summary: Nullable<string>
  gpt_summary: Nullable<string>
  created: CreatedTimestamp
}

export interface FPLTable {
  id: string
  process: string
  parent: Nullable<string>
  playbook_metadata: Nullable<string>
  cell_metadata: Nullable<string>
  created: CreatedTimestamp
  hits: Counter
}

// user content

export interface UserPlaybookTable {
  id: string
  user: string
  playbook: string
  title: Nullable<string>
  description: Nullable<string>
  inputs: Nullable<string>
  outputs: Nullable<string>
  bco: Nullable<string>
  clicks: Counter
  public: ColumnType<boolean, boolean | undefined, boolean>
  created: CreatedTimestamp
}

export interface UploadTable {
  id: string
  url: string
  sha256: string
  size: number
  created: CreatedTimestamp
}

export interface UserUploadTable {
  id: string
  user: Nullable<string>
  upload: string
  filename: string
  created: CreatedTimestamp
}

export interface ThreadTable {
  id: string
  user: Nullable<string>
  openai_thread: Nullable<string>
  created: CreatedTimestamp
}

export interface ThreadMessageTable {
  id: string
  fpl: Nullable<string>
  thread: string
  role: string
  content: string
  feedback: Nullable<string>
  created: CreatedTimestamp
}

export interface Database {
  user: UserTable
  account: AccountTable
  session: SessionTable
  verification_token: VerificationTokenTable
  user_integrations: UserIntegrationsTable
  proxy_session: ProxySessionTable
  suggestion: SuggestionTable
  data: DataTable
  process: ProcessTable
  process_input: ProcessInputTable
  resolved: ResolvedTable
  cell_metadata: CellMetadataTable
  playbook_metadata: PlaybookMetadataTable
  fpl: FPLTable
  user_playbook: UserPlaybookTable
  upload: UploadTable
  user_upload: UserUploadTable
  thread: ThreadTable
  thread_message: ThreadMessageTable
}

export type User = Selectable<UserTable>
export type NewUser = Insertable<UserTable>
export type UserUpdate = Updateable<UserTable>
export type Account = Selectable<AccountTable>
export type Session = Selectable<SessionTable>
export type VerificationToken = Selectable<VerificationTokenTable>
export type UserIntegrations = Selectable<UserIntegrationsTable>
export type ProxySession = Selectable<ProxySessionTable>
export type Suggestion = Selectable<SuggestionTable>
export type UserPlaybook = Selectable<UserPlaybookTable>
export type Upload = Selectable<UploadTable>
export type UserUpload = Selectable<UserUploadTable>
export type Thread = Selectable<ThreadTable>
export type ThreadMessage = Selectable<ThreadMessageTable>
/** user_upload joined with its upload */
export type UserUploadComplete = Pick<UserUpload, 'id' | 'user' | 'filename' | 'created'> & Pick<Upload, 'url' | 'sha256' | 'size'>
