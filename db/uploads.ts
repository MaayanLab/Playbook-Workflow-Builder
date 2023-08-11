import { Table, View } from '@/spec/sql'
import { z } from 'zod'
import { v4 as uuidv4 } from 'uuid'
import { z_uuid, z_bigint_codec } from '@/utils/zod'

export const upload = Table.create('upload')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('url', 'varchar', 'not null', z.string())
  .field('sha256', 'varchar', 'not null', z.string())
  .field('size', 'bigint', 'not null', z_bigint_codec())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
  .build()

export const user_upload = Table.create('user_upload')
  .field('id', 'uuid', 'default uuid_generate_v4()', z_uuid(), { primaryKey: true, default: uuidv4 })
  .field('user', 'uuid', 'references "user" ("id") on delete cascade', z_uuid().nullable())
  .field('upload', 'uuid', 'not null references "upload" ("id") on delete cascade', z_uuid())
  .field('filename', 'varchar', 'not null', z.string())
  .field('created', 'timestamp', 'not null default now()', z.date(), { default: () => new Date() })
  .build()

export const user_upload_complete = View.create('user_upload_complete')
  .field('id', z_uuid(), { primaryKey: true })
  .field('user', z_uuid())
  .field('url', z.string())
  .field('sha256', z.string())
  .field('size', z_bigint_codec())
  .field('filename', z.string())
  .field('created', z.date())
  .sql(`
    select
      "user_upload"."id",
      "user_upload"."user",
      "upload"."url",
      "upload"."sha256",
      "upload"."size",
      "user_upload"."filename",
      "user_upload"."created"
    from "user_upload"
    left join "upload" on "user_upload"."upload" = "upload"."id";
  `)
  .js(async (db: any) => await Promise.all((await db.objects.user_upload.findMany()).map(async (user_upload: any) => {
    const upload = await db.objects.upload.findUnique({ where: { id: user_upload.upload } })
    return {
      id: user_upload.id,
      user: user_upload.user,
      url: upload.url,
      sha256: upload.sha256,
      size: upload.size,
      filename: user_upload.filename,
      created: user_upload.created,
    }
  })))
  .build()
