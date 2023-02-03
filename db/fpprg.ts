import { SQL, Table, View } from '@/spec/sql'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import { z_uuid } from '@/utils/zod'

export const data = Table.create('data')
  .field('id', 'uuid', '', z_uuid(), { primaryKey: true })
  .field('type', 'varchar', 'not null', z.string())
  .field('value', 'varchar', 'not null', z.string())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
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
  .field('id', 'uuid', '', z_uuid(), { primaryKey: true })
  .field('type', 'varchar', 'not null', z.string())
  .field('data', 'uuid', 'references data ("id") on delete cascade', z_uuid().nullable())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
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
  .field('id', 'uuid', 'references process ("id") on delete cascade', z_uuid(), { primaryKey: true })
  .field('key', 'varchar', 'not null', z.string(), { primaryKey: true })
  .field('value', 'uuid', 'not null references process ("id") on delete cascade', z_uuid())
  .build()

export const resolved = Table.create('resolved')
  .field('id', 'uuid', 'references process ("id") on delete cascade', z_uuid(), { primaryKey: true })
  .field('data', 'uuid', 'references data ("id") on delete cascade', z_uuid().nullable())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
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

export const process_complete = View.create('process_complete')
    .field('id', z_uuid(), { primaryKey: true })
    .field('type', z.string())
    .field('data', z.string().nullable())
    .field('inputs', z.record(z.string(), z_uuid()))
    .field('resolved', z.boolean())
    .field('output', z_uuid().nullable())
    .sql(`
      select
        "process"."id",
        "process"."type",
        "process"."data",
        coalesce(
          (select jsonb_object_agg("process_input"."key", "process_input"."value")
           from "process_input"
           where process."id" = "process_input"."id"),
          '{}'::jsonb
        ) as "inputs",
        ("resolved"."id" is not null) as "resolved",
        "resolved"."data" as "output"
      from "process"
      left join "resolved" on "process"."id" = "resolved"."id";
    `)
    .js(async (db: any) => await Promise.all((await db.objects.process.findMany()).map(async ({ id, type, data }: any) => {
      const resolved = await db.objects.resolved.findUnique({ where: { id } })
      const inputs = dict.init(await db.objects.process_input.findMany({ where: { id } }))
      return {
        id,
        type,
        data,
        inputs,
        resolved: resolved !== null,
        output: resolved !== null ? resolved.data : null,
      }
    })))
    .build()

export const fpl = Table.create('fpl')
  .field('id', 'uuid', '', z_uuid(), { primaryKey: true })
  .field('process', 'uuid', 'not null references process ("id") on delete cascade', z_uuid())
  .field('parent', 'uuid', 'references fpl ("id") on delete cascade', z_uuid().nullable())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
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
