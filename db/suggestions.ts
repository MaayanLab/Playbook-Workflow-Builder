import { Table } from '@/spec/sql'
import * as z from 'zod'

const z_uuid = z.string

export const suggestion = Table.create('suggestion')
  .field('id', 'uuid', 'primary key default uuid_generate_v4()', z_uuid())
  .field('name', 'varchar', 'not null', z.string())
  .field('inputs', 'varchar', 'not null', z.string())
  .field('output', 'varchar', 'not null', z.string())
  .field('author_name', 'varchar', 'not null', z.string())
  .field('author_email', 'varchar', 'not null', z.string())
  .field('author_org', 'varchar', 'not null', z.string())
  .field('description', 'varchar', 'not null', z.string())
  .build()
