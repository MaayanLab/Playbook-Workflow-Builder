import { Table } from '@/spec/sql'
import * as z from 'zod'

const z_uuid = z.string

export const user = Table.create('user')
  .field('id', 'uuid', 'primary key', z_uuid())
  .field('name', 'varchar', '', z.string().optional())
  .field('email', 'varchar', '', z.string().optional())
  .field('emailVerified', 'timestamp', '', z.date())
  .field('image', 'varchar', '', z.string().optional())
  .build()


export const account = Table.create('account')
  .field('id', 'uuid', 'primary key', z_uuid())
  .field('userId', 'uuid', 'not null references user ("id") on delete cascade', z_uuid())
  .field('type', 'varchar', '', z.string().optional())
  .field('provider', 'varchar', '', z.string().optional())
  .field('providerAccountId', 'varchar', '', z.string().optional())
  .field('refresh_token', 'varchar', '', z.string().optional())
  .field('access_token', 'varchar', '', z.string().optional())
  .field('expires_at', 'int', '', z.number())
  .field('token_type', 'varchar', '', z.string().optional())
  .field('scope', 'varchar', '', z.string().optional())
  .field('id_token', 'varchar', '', z.string().optional())
  .field('session_state', 'varchar', '', z.string().optional())
  .field('oauth_token_secret', 'varchar', '', z.string().optional())
  .field('oauth_token', 'varchar', '', z.string().optional())
  .build()

export const session = Table.create('session')
  .field('id', 'uuid', 'primary_key', z_uuid())
  .field('expires', 'timestamp', '', z.date())
  .field('sessionToken', 'varchar', '', z.string().optional())
  .field('userId', 'uuid', 'not null references user ("id") on delete cascade', z_uuid())
  .build()

export const verification_token = Table.create('verification_token')
  .field('identifier', 'uuid', 'primary key', z_uuid())
  .field('token', 'varchar', '', z.string().optional())
  .field('expires', 'timestamp', '', z.date())
  .build()
