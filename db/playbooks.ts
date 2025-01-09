import { Table } from '@/spec/sql'
import { z } from 'zod'
import { v4 as uuidv4 } from 'uuid'
import { z_uuid, nullable_undefined_codec, z_bigint_codec, json_safe_timestamp_codec } from '@/utils/zod'

export const user_playbook = Table.create('user_playbook')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('user', 'uuid', 'not null references "user" ("id") on delete cascade', z_uuid())
  .field('playbook', 'uuid', 'not null references "fpl" ("id") on delete cascade', z_uuid())
  .field('title', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('description', 'text', '', nullable_undefined_codec(z.string()))
  .field('inputs', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('outputs', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('bco', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('clicks', 'bigint', 'not null default 0', z_bigint_codec(), { default: () => 0 })
  .field('public', 'boolean', 'not null default false', z.boolean(), { default: () => false })
  .field('created', 'timestamp without time zone', 'not null default now()', json_safe_timestamp_codec(), { default: () => new Date() })
  .build()
