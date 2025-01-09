import { Table } from '@/spec/sql'
import { z } from 'zod'
import { v4 as uuidv4 } from 'uuid'
import { json_safe_timestamp_codec, z_uuid } from '@/utils/zod'

export const proxy_session = Table.create('proxy_session')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('state', 'varchar', '', z.string().nullable(), { default: () => null })
  .field('run_id', 'varchar', '', z.string().nullable(), { default: () => null })
  .field('created', 'timestamp without time zone', 'not null default now()', json_safe_timestamp_codec(), { default: () => new Date() })
  .build()
