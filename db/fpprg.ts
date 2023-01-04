import { Table, View } from '@/spec/sql'
import * as z from 'zod'

const z_uuid = z.string

export const data = Table.create('data')
  .field('id', 'uuid', 'primary key', z_uuid())
  .field('type', 'varchar', '', z.string())
  .field('value', 'jsonb', '', z.record(z.any(), z.any()).optional())
  .extra_insert('on conflict ("id") do nothing')
  .build()

export const process = Table.create('process')
  .field('id', 'uuid', 'primary key', z_uuid())
  .field('type', 'varchar', 'not null', z.string())
  .field('data', 'uuid', 'foreign key references data ("id") on delete cascade', z_uuid().optional())
  .extra_insert('on conflict ("id") do nothing')
  .build()

export const process_input = Table.create('process_input')
  .field('id', 'uuid', 'foreign key references process ("id") on delete cascade', z_uuid())
  .field('key', 'varchar', 'not null', z.string())
  .field('value', 'uuid', 'not null foreign key references process ("id") on delete cascade', z_uuid())
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
      jsonb_object_agg(process_inputs."key", process_inputs."value") as "inputs"
    from process
    left join process_inputs on process."id" = process_inputs."id"
    group by process."id";
  `)
  .build()

export const resolved = Table.create('resolved')
  .field('id', 'uuid', 'primary key foreign key references process ("id") on delete cascade', z_uuid())
  .field('data', 'uuid', 'foreign key references data ("id") on delete cascade', z_uuid().optional())
  .extra_insert('on conflict ("id") do nothing')
  .build()

export const fpl = Table.create('fpl')
  .field('id', 'uuid', 'primary key', z_uuid())
  .field('process', 'uuid', 'not null foreign key references process ("id") on delete cascade', z_uuid())
  .field('parent', 'uuid', 'foreign key references fpl ("id") on delete cascade', z_uuid().optional())
  .extra_insert('on conflict ("id") do nothing')
  .build()
