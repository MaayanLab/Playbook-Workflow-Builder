import { Table } from '@/spec/sql'
import { z } from 'zod'
import { z_uuid } from '@/utils/zod'

export const fpl_embedding = Table.create('fpl_embedding')
  .field('id', 'uuid', 'primary key references "fpl" ("id") on delete cascade', z_uuid())
  .field('embedding', 'vector(1536)', 'not null', z.number().array())
  .build()
