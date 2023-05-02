import { SQL } from '@/spec/sql'
import { z } from 'zod'

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
          'id', rec.id
        )::text
      );

      return rec;
    end; $$ language plpgsql;
  `)
  .down(`
    drop function notify_trigger() cascade;
  `)
  .build()

const z_notify_insertion_trigger_payload = z.object({
  'timestamp': z.string(),
  'operation': z.string(),
  'schema': z.string(),
  'table': z.string(),
  'label': z.array(z.string()),
  'id': z.string(),
})

export const notify_insertion_trigger_payload = {
  codec: {
    encode: z_notify_insertion_trigger_payload.parse,
    decode: z_notify_insertion_trigger_payload.parse,
  }
}

export const assert_is_one = SQL.create()
  .up(`
    create function assert_is_one(value bigint) returns void as $$
    begin
      if value != 1 then
        raise 'Expected one got %', value using errcode = 'unique_violation';
      end if;
    end;
    $$ language plpgsql;
  `)
  .down(`drop function assert_is_one(bigint);`)
  .build()
