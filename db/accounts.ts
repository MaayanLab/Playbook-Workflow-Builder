import { SQL, Table } from '@/spec/sql'
import { z } from 'zod'
import { v4 as uuidv4 } from 'uuid'
import { z_uuid, nullable_undefined_codec, z_bigint_codec } from '@/utils/zod'


export const user = Table.create('user')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('name', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('email', 'varchar', 'not null', z.string())
  .field('emailVerified', 'timestamp', '', z.date().nullable())
  .field('image', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('affiliation', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
  .build()

export const user_email_index = SQL.create()
  .up('create index user_email_idx on "user" ("email");')
  .down('drop index user_email_idx;')
  .build()

export const account = Table.create('account')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('userId', 'uuid', 'not null references "user" ("id") on delete cascade', z_uuid())
  .field('type', 'varchar', 'not null', z.enum(["oauth", "email", "credentials"]))
  .field('provider', 'varchar', 'not null', z.string())
  .field('providerAccountId', 'varchar', 'not null', z.string())
  .field('refresh_token', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('access_token', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('expires_at', 'bigint', '', z_bigint_codec())
  .field('token_type', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('scope', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('id_token', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('session_state', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('oauth_token_secret', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('oauth_token', 'varchar', '', nullable_undefined_codec(z.string()))
  .build()

export const session = Table.create('session')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('expires', 'timestamp', '', z.date())
  .field('sessionToken', 'varchar', 'not null', z.string())
  .field('userId', 'uuid', 'not null references "user" ("id") on delete cascade', z_uuid())
  .build()

export const verification_token = Table.create('verification_token')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('identifier', 'varchar', '', z.string())
  .field('token', 'varchar', '', z.string())
  .field('expires', 'timestamp', '', z.date())
  .build()

export const user_integrations = Table.create('user_integrations')
  .field('id', 'uuid', 'not null references "user" ("id") on delete cascade', z_uuid(), { primaryKey: true })
  .field('cavatica_api_key', 'varchar', '', z.string()) 
  .field('cavatica_default_project', 'varchar', '', z.string())
  .build()
