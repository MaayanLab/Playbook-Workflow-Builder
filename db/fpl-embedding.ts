import { SQL, Table } from '@/spec/sql'
import { z } from 'zod'
import { z_uuid } from '@/utils/zod'

export const pg_vector = SQL.create()
  .up(`create extension if not exists "vector";`)
  .down(`drop extension "vector";`)
  .build()

export const fpl_embedding = Table.create('fpl_embedding')
  .field('id', 'uuid', 'references "fpl" ("id") on delete cascade', z_uuid(), { primaryKey: true })
  .field('embedding', 'vector(1536)', 'not null', z.number().array())
  .build()
