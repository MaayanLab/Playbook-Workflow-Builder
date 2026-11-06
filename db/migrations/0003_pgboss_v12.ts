/**
 * pg-boss 9 → 12
 *
 * pg-boss manages its own schema (`pgboss`) and migrates it on start, but it can't migrate
 *  a pg-boss 9 schema (schema version 20): pg-boss 10 dropped support for migrating from 9,
 *  and pg-boss 12 can only migrate from schema version 25 (pg-boss 11).
 *
 * The old schema is renamed to `pgboss_v9` so that pg-boss 12 creates a fresh `pgboss` schema
 *  when it starts. Nothing of value is lost: the queue only holds requests to resolve
 *  processes, which are re-submitted while a client awaits the result.
 *  `pgboss_v9` (job history) can be dropped once the upgrade is confirmed.
 *
 * Stop all ui & worker processes running pg-boss 9 before applying this migration.
 */
import { type Kysely, PostgresAdapter, sql } from 'kysely'

export async function up(db: Kysely<any>): Promise<void> {
  if (!(db.getExecutor().adapter instanceof PostgresAdapter)) return
  await sql`
    do $$
    declare
      schema_version int;
    begin
      if exists (select 1 from information_schema.tables where table_schema = 'pgboss' and table_name = 'version') then
        execute 'select max(version) from pgboss.version' into schema_version;
        if schema_version is null or schema_version < 25 then
          alter schema pgboss rename to pgboss_v9;
        end if;
      end if;
    end $$;
  `.execute(db)
}

export async function down(db: Kysely<any>): Promise<void> {
  if (!(db.getExecutor().adapter instanceof PostgresAdapter)) return
  await sql`
    do $$
    begin
      if exists (select 1 from information_schema.schemata where schema_name = 'pgboss_v9') then
        drop schema if exists pgboss cascade;
        alter schema pgboss_v9 rename to pgboss;
      end if;
    end $$;
  `.execute(db)
}
