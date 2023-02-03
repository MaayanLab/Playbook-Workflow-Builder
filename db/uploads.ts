import { Table } from '@/spec/sql'
import { z } from 'zod'
import { v4 as uuidv4 } from 'uuid'
import { z_uuid } from '@/utils/zod'

export const upload = Table.create('upload')
  .field('url', 'varchar', 'not null', z.string(), { primaryKey: true })
  .field('sha256', 'varchar', 'not null', z.string())
  .field('size', 'number', 'not null', z.number())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
  .build()

export const user_upload = Table.create('user_upload')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('user', 'uuid', 'not null references "user" ("id") on delete cascade', z_uuid())
  .field('upload', 'uuid', 'not null references "upload" ("id") on delete cascade', z_uuid())
  .field('filename', 'varchar', 'not null', z.string())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
  .build()
