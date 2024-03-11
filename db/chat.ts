import { Table } from '@/spec/sql'
import { z } from 'zod'
import { v4 as uuidv4 } from 'uuid'
import { z_uuid } from '@/utils/zod'

export const chat_message = Table.create('chat_message')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('user', 'uuid', 'not null references "user" ("id") on delete cascade', z_uuid())
  .field('previous', 'uuid', '', z_uuid())
  .field('role', 'varchar', 'not null', z.string())
  .field('content', 'text', 'not null', z.string())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
  .build()
