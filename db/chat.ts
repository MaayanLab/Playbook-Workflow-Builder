import { Table } from '@/spec/sql'
import { z } from 'zod'
import { v4 as uuidv4 } from 'uuid'
import { z_uuid } from '@/utils/zod'

export const thread = Table.create('thread')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('user', 'uuid', 'references "user" ("id") on delete cascade', z_uuid())
  .field('openai_thread', 'varchar', 'not null', z.string())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
  .build()

export const thread_message = Table.create('thread_message')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('fpl', 'uuid', 'references "fpl" ("id") on delete cascade', z_uuid().nullable())
  .field('thread', 'uuid', 'not null references "thread" ("id") on delete cascade', z_uuid())
  .field('role', 'varchar', 'not null', z.string())
  .field('content', 'varchar', 'not null', z.string())
  .field('feedback', 'varchar', '', z.string().nullable())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
  .build()
