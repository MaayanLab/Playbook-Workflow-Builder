import { SQL, Table } from '@/spec/sql'
import * as z from 'zod'

const z_uuid = z.string
const nullable_undefined_codec = <C>(type: z.ZodType<C>) => ({
  decode: type.nullable().transform(v => v !== null ? v : undefined).parse,
  encode: type.optional().transform(v => v !== undefined ? v : null).parse,
})

export const user = Table.create('user')
  .field('id', 'uuid', 'primary key default uuid_generate_v4()', z_uuid())
  .field('name', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('email', 'varchar', 'not null', z.string())
  .field('emailVerified', 'timestamp', '', z.date())
  .field('image', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('affiliation', 'varchar', '', nullable_undefined_codec(z.string()))
  .build()

export const user_email_index = SQL.create()
  .up('create index user_email_idx on "user" ("email");')
  .down('drop index user_email_idx;')
  .build()

export const account = Table.create('account')
  .field('id', 'uuid', 'primary key default uuid_generate_v4()', z_uuid())
  .field('userId', 'uuid', 'not null references "user" ("id") on delete cascade', z_uuid())
  .field('type', 'varchar', 'not null', z.enum(["oauth", "email", "credentials"]))
  .field('provider', 'varchar', 'not null', z.string())
  .field('providerAccountId', 'varchar', 'not null', z.string())
  .field('refresh_token', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('access_token', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('expires_at', 'int', '', z.number())
  .field('token_type', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('scope', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('id_token', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('session_state', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('oauth_token_secret', 'varchar', '', nullable_undefined_codec(z.string()))
  .field('oauth_token', 'varchar', '', nullable_undefined_codec(z.string()))
  .build()

export const session = Table.create('session')
  .field('id', 'uuid', 'primary key default uuid_generate_v4()', z_uuid())
  .field('expires', 'timestamp', '', z.date())
  .field('sessionToken', 'varchar', 'not null', z.string())
  .field('userId', 'uuid', 'not null references "user" ("id") on delete cascade', z_uuid())
  .build()

export const verification_token = Table.create('verification_token')
  .field('identifier', 'varchar', 'primary key', z.string())
  .field('token', 'varchar', '', z.string())
  .field('expires', 'timestamp', '', z.date())
  .build()
