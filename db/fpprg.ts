import { SQL, Table, View } from '@/spec/sql'
import * as z from 'zod'

const z_uuid = z.string

export const pg_uuids = SQL.create()
  .up(`create extension if not exists "uuid-ossp";`)
  .down(`drop extension "uuid-ossp";`)
  .build()

export const notify_insertion_trigger = SQL.create()
  .up(`
    create or replace function notify_trigger() returns trigger as $$
    declare
      rec record;
    begin
      case TG_OP
        when 'INSERT', 'UPDATE' then rec := NEW;
        when 'DELETE' then rec := OLD;
        else raise exception 'Unhandled TG_OP: "%"', TG_OP;
      end case;

      perform pg_notify(
        'on_insert',
        json_build_object(
          'timestamp', CURRENT_TIMESTAMP,
          'operation', TG_OP,
          'schema', TG_TABLE_SCHEMA,
          'table', TG_TABLE_NAME,
          'label', array_to_json(TG_ARGV),
          'body', row_to_json(rec)
        )::text
      );

      return rec;
    end; $$ language plpgsql;
  `)
  .down(`
    drop function notify_trigger() cascade;
  `)
  .build()

export const notify_insertion_trigger_payload = z.object({
  'timestamp': z.string(),
  'operation': z.string(),
  'schema': z.string(),
  'table': z.string(),
  'label': z.array(z.string()),
  'body': z.record(z.string(), z.any()),
})

export const data = Table.create('data')
  .field('id', 'uuid', 'primary key', z_uuid())
  .field('type', 'varchar', 'not null', z.string())
  .field('value', 'varchar', 'not null', z.string())
  .extra_insert('on conflict ("id") do nothing')
  .build()

export const data_trigger = SQL.create()
  .up(`
    create trigger data_notify after insert on data
    for each row execute procedure notify_trigger('data');
  `)
  .down(`
    drop trigger data_notify on data cascade;
  `)
  .build()

export const process = Table.create('process')
  .field('id', 'uuid', 'primary key', z_uuid())
  .field('type', 'varchar', 'not null', z.string())
  .field('data', 'uuid', 'references data ("id") on delete cascade', z_uuid().optional())
  .extra_insert('on conflict ("id") do nothing')
  .build()

export const process_trigger = SQL.create()
  .up(`
    create trigger process_notify after insert on process
    for each row execute procedure notify_trigger('process');
  `)
  .down(`
    drop trigger process_notify on process cascade;
  `)
  .build()

export const process_input = Table.create('process_input')
  .field('id', 'uuid', 'references process ("id") on delete cascade', z_uuid())
  .field('key', 'varchar', 'not null', z.string())
  .field('value', 'uuid', 'not null references process ("id") on delete cascade', z_uuid())
  .extra('primary key ("id", "key")')
  .extra_insert('on conflict ("id", "key") do nothing')
  .build()

export const process_complete = View.create('process_complete')
  .field('id', z_uuid())
  .field('type', z.string())
  .field('data', z.string().optional())
  .field('inputs', z.record(z.string(), z_uuid()))
  .sql(`
    select
      process."id",
      process."type",
      process."data",
      jsonb_object_agg(process_input."key", process_input."value") as "inputs"
    from process
    left join process_input on process."id" = process_input."id"
    group by process."id";
  `)
  .build()

export const resolved = Table.create('resolved')
  .field('id', 'uuid', 'primary key references process ("id") on delete cascade', z_uuid())
  .field('data', 'uuid', 'references data ("id") on delete cascade', z_uuid().optional())
  .extra_insert('on conflict ("id") do nothing')
  .build()

export const resolved_trigger = SQL.create()
  .up(`
    create trigger resolved_notify after insert on resolved
    for each row execute procedure notify_trigger('resolved');
  `)
  .down(`
    drop trigger resolved_notify on resolved cascade;
  `)
  .build()

export const fpl = Table.create('fpl')
  .field('id', 'uuid', 'primary key', z_uuid())
  .field('process', 'uuid', 'not null references process ("id") on delete cascade', z_uuid())
  .field('parent', 'uuid', 'references fpl ("id") on delete cascade', z_uuid().optional())
  .extra_insert('on conflict ("id") do nothing')
  .build()

export const fpl_trigger = SQL.create()
  .up(`
    create trigger fpl_notify after insert on fpl
    for each row execute procedure notify_trigger('fpl');
  `)
  .down(`
    drop trigger fpl_notify on fpl cascade;
  `)
  .build()
