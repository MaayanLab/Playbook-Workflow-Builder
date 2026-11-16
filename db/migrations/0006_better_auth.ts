/** Better Auth's Kysely schema, reusing the application's existing user table. */
import { type Kysely, PostgresAdapter, sql } from 'kysely'

export async function up(db: Kysely<any>): Promise<void> {
  const pg = db.getExecutor().adapter instanceof PostgresAdapter
  const text = 'text' as const
  const timestamp = pg ? 'timestamp' : 'text'
  const boolean = 'boolean' as const
  const now = pg ? sql`now()` : sql`(strftime('%Y-%m-%dT%H:%M:%fZ', 'now'))`

  await db.schema.alterTable('user')
    .addColumn('authEmailVerified', boolean, col => col.notNull().defaultTo(pg ? false : 0))
    .execute()
  await db.schema.alterTable('user')
    .addColumn('authUpdated', timestamp)
    .execute()
  await db.schema.alterTable('user')
    .addColumn('isAnonymous', boolean, col => col.notNull().defaultTo(pg ? false : 0))
    .execute()
  await sql`update "user" set "authUpdated" = ${now} where "authUpdated" is null`.execute(db)
  await sql`update "user" set "authEmailVerified" = true where "emailVerified" is not null`.execute(db)

  await db.schema.createTable('auth_session')
    .addColumn('id', text, col => col.primaryKey())
    .addColumn('token', text, col => col.notNull().unique())
    .addColumn('expiresAt', timestamp, col => col.notNull())
    .addColumn('createdAt', timestamp, col => col.notNull().defaultTo(now))
    .addColumn('updatedAt', timestamp, col => col.notNull().defaultTo(now))
    .addColumn('ipAddress', text)
    .addColumn('userAgent', text)
    .addColumn('userId', pg ? 'uuid' : text, col => col.notNull().references('user.id').onDelete('cascade'))
    .execute()

  await db.schema.createTable('auth_account')
    .addColumn('id', text, col => col.primaryKey())
    .addColumn('accountId', text, col => col.notNull())
    .addColumn('providerId', text, col => col.notNull())
    .addColumn('userId', pg ? 'uuid' : text, col => col.notNull().references('user.id').onDelete('cascade'))
    .addColumn('accessToken', text)
    .addColumn('refreshToken', text)
    .addColumn('idToken', text)
    .addColumn('accessTokenExpiresAt', timestamp)
    .addColumn('refreshTokenExpiresAt', timestamp)
    .addColumn('scope', text)
    .addColumn('password', text)
    .addColumn('createdAt', timestamp, col => col.notNull().defaultTo(now))
    .addColumn('updatedAt', timestamp, col => col.notNull().defaultTo(now))
    .execute()
  await db.schema.createIndex('auth_account_provider_idx').on('auth_account').columns(['providerId', 'accountId']).unique().execute()
  await sql`
    insert into "auth_account" (
      "id", "accountId", "providerId", "userId", "accessToken", "refreshToken", "idToken", "scope", "createdAt", "updatedAt"
    )
    select "id", "providerAccountId", "provider", "userId", "access_token", "refresh_token", "id_token", "scope", ${now}, ${now}
    from "account"
    where true
    on conflict ("providerId", "accountId") do nothing
  `.execute(db)
  if (pg) {
    await sql`
      update "auth_account" as target
      set "accessTokenExpiresAt" = to_timestamp(source."expires_at")
      from "account" as source
      where target."id" = source."id"::text and source."expires_at" is not null
    `.execute(db)
  } else {
    await sql`
      update "auth_account"
      set "accessTokenExpiresAt" = (
        select strftime('%Y-%m-%dT%H:%M:%fZ', source."expires_at", 'unixepoch')
        from "account" as source
        where source."id" = "auth_account"."id"
      )
      where exists (
        select 1 from "account" as source
        where source."id" = "auth_account"."id" and source."expires_at" is not null
      )
    `.execute(db)
  }

  await db.schema.createTable('auth_verification')
    .addColumn('id', text, col => col.primaryKey())
    .addColumn('identifier', text, col => col.notNull())
    .addColumn('value', text, col => col.notNull())
    .addColumn('expiresAt', timestamp, col => col.notNull())
    .addColumn('createdAt', timestamp, col => col.notNull().defaultTo(now))
    .addColumn('updatedAt', timestamp, col => col.notNull().defaultTo(now))
    .execute()
}

export async function down(db: Kysely<any>): Promise<void> {
  await db.schema.dropTable('auth_verification').execute()
  await db.schema.dropTable('auth_account').execute()
  await db.schema.dropTable('auth_session').execute()
  await db.schema.alterTable('user').dropColumn('authUpdated').execute()
  await db.schema.alterTable('user').dropColumn('isAnonymous').execute()
  await db.schema.alterTable('user').dropColumn('authEmailVerified').execute()
}
