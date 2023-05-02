import { Table } from '@/spec/sql'
import { z } from 'zod'
import { v4 as uuidv4 } from 'uuid'
import { json_safe_timestamp_codec, z_uuid } from '@/utils/zod'

export const suggestion = Table.create('suggestion')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('user', 'uuid', 'not null references "user" ("id") on delete cascade', z_uuid())
  .field('name', 'varchar', 'not null', z.string())
  .field('inputs', 'varchar', 'not null', z.string())
  .field('output', 'varchar', 'not null', z.string())
  .field('description', 'varchar', 'not null', z.string())
  .field('created', 'timestamp', 'not null default now()', json_safe_timestamp_codec(), { default: () => (new Date()).toISOString() })
  .build()
