/**
 * A lightweight "ORM" for constructing type-safe pg tables, useful so
 *  we can define the .sql and the typescript model at the same time.
 **/
import { z } from "zod"
import * as dict from '@/utils/dict'
import { Codec, Decoded, Encoded } from "./codec"

const identityZodCodec = <T>(z: z.ZodType<T>): Codec<T, T> => ({ encode: z.parse, decode: z.parse })

const ObjectCodec = <T>(field_codecs: {[K in keyof T]: Codec<Decoded<T[K]>, Encoded<T[K]>>}) => ({
  encode: (data: { [K in keyof T]: Decoded<T[K]> }) =>
    dict.init(dict.items(data).map(({ key, value }) => ({ key, value: field_codecs[key].encode(value) }))),
  decode: (data: { [K in keyof T]: Encoded<T[K]> }) =>
    dict.init(dict.items(data).map(({ key, value }) => ({ key, value: field_codecs[key].decode(value) }))),
})

export type TypedSchema<T = {}> = {
  name: string
  codec: Codec<{[K in keyof T]: Decoded<T[K]>}, {[K in keyof T]: Encoded<T[K]>}>
  field_pk: Array<keyof T>
  field_default: Partial<{[K in keyof T]: Decoded<T[K]>}>
  field_codecs: {[K in keyof T]: Codec<Decoded<T[K]>, Encoded<T[K]>>}
  field_sql: {[K in keyof T]: string}
  schema_up: string
  schema_down: string
  js?: (db: any) => Promise<Array<{[K in keyof T]: Decoded<T[K]>}>>
}
export type TypedSchemaRecord<Schema> = Schema extends TypedSchema<infer T> ? {[K in keyof T]: Decoded<T[K]>} : never

export type TableSchema<T = {}> = {
  name: string
  field_pk: {[K in keyof T]: boolean}
  field_default: {[K in keyof T]?: () => Decoded<T[K]>}
  field_codecs: {[K in keyof T]: Codec<Decoded<T[K]>, Encoded<T[K]>>}
  field_sql: {[K in keyof T]: string}
  field_extra_sql: {[K in keyof T]: string}
  extra_sql: string[]
}

/**
 * Usage:
 * const table_name = Table.createTable('table_name')
 *   // fieldname, pg schema, and a codec to determine how it should be treated in and out of the db
 *   .field('fieldname', 'varchar', { encode: z.string().transform(s => s|0).parse, decode: z.int().transform(i => i+'').parse })
 *   // in the case of storing it equally both ways, a zod schema can be used
 *   .field('fieldname', 'varchar', z.string())
 *   .extra('') // any extra sql in the table definition like primary key/constraints
 *   .build()
 */
export class Table<T = {}> {
  constructor(public t: TableSchema<T>) {}
  static create(name: string) {
    return new Table({
      name,
      field_codecs: {},
      field_pk: {},
      field_default: {},
      field_sql: {},
      field_extra_sql: {},
      extra_sql: [],
    })
  }
  field<S extends string, D, E = D>(name: S, sql: string, extra_sql: string, type: z.ZodType<D> | Codec<D, E>, opts?: { primaryKey?: boolean, default?: () => D }) {
    return new Table({
      ...this.t,
      field_pk: {...this.t.field_pk, [name]: !!(opts||{}).primaryKey },
      field_default: {...this.t.field_default, [name]: (opts||{}).default },
      field_codecs: { ...this.t.field_codecs, [name]: 'parse' in type ? identityZodCodec(type): type },
      field_sql: { ...this.t.field_sql, [name]: sql },
      field_extra_sql: { ...this.t.field_extra_sql, [name]: extra_sql },
    } as TableSchema<T & { [K in S]: Codec<D, E> }>)
  }
  extra(extra_sql: string) {
    return new Table({
      ...this.t,
      extra_sql: [...this.t.extra_sql, extra_sql],
    })
  }
  build() {
    const { name, field_codecs, field_sql, field_extra_sql } = this.t
    const field_pk = dict.items(this.t.field_pk).filter(({ value }) => !!value).map(({ key }) => key)
    if (field_pk.length <= 0) throw new Error(`Missing primary key on ${name}`)
    const field_default = dict.init(dict.items(this.t.field_default).filter(({ value }) => !!value))
    const codec = ObjectCodec(field_codecs)
    const schema_up = [
      `create table ${JSON.stringify(this.t.name)} (`,
      [
        ...dict.sortedKeys(field_sql).map(field =>
          `  ${JSON.stringify(field)} ${field_sql[field]} ${field_extra_sql[field]}`
        ),
        `primary key (${field_pk.map((key) => JSON.stringify(key)).join(', ')})`,
        ...this.t.extra_sql.map(sql => `  ${sql}`),
      ].join(',\n'),
      `);`,
    ].join('\n')

    const schema_down = `drop table ${JSON.stringify(this.t.name)};`
    return { codec, name, field_sql, field_pk, field_default, field_codecs, schema_up, schema_down } as TypedSchema<T>
  }
}

export type ViewSchema<T = {}> = {
  name: string
  field_pk: {[K in keyof T]: boolean}
  field_codecs: {[K in keyof T]: Codec<Decoded<T[K]>, Encoded<T[K]>>}
  sql: string
  js: (db: any) => any
}

export class View<T extends {} = {}> {
  constructor(public t: ViewSchema<T>) {}
  static create(name: string) {
    return new View({ name, field_codecs: {}, field_pk: {}, sql: '', js: () => [] })
  }
  field<S extends string, D, E = D>(name: S, type: z.ZodType<D> | Codec<D, E>, opts?: { primaryKey?: boolean, default?: () => D }) {
    return new View({
      ...this.t,
      field_pk: {...this.t.field_pk, [name]: !!(opts||{}).primaryKey },
      field_codecs: { ...this.t.field_codecs, [name]: 'parse' in type ? identityZodCodec(type): type },
    } as ViewSchema<T & { [K in S]: Codec<D, E> }>)
  }
  /**
   * SQL Version of the view
   */
  sql(...sql: string[]) {
    return new View({ ...this.t, sql: sql.join('\n'), })
  }
  /**
   * Construct the view in-memory with JS (not particularly efficient but it's not going to be used in prod)
   */
  js(js: (db: any) => Promise<any>) {
    return new View({ ...this.t, js, })
  }
  build() {
    const { name, field_codecs, js } = this.t
    const field_pk = dict.items(this.t.field_pk).filter(({ value }) => !!value).map(({ key }) => key)
    if (field_pk.length <= 0) throw new Error(`Missing primary key on ${name}`)
    const codec = ObjectCodec(field_codecs)
    const schema_up = `create view ${JSON.stringify(this.t.name)} as ${this.t.sql}`
    const schema_down = `drop view ${JSON.stringify(this.t.name)};`
    return { codec, name, js, field_pk, field_codecs, schema_up, schema_down } as TypedSchema<T>
  }
}


export type SQLSchema = {
  schema_up: string
  schema_down: string
}

export class SQL {
  constructor(public t: SQLSchema) {}
  static create() {
    return new SQL({ schema_up: '', schema_down: '' })
  }
  up(...sql: string[]) {
    return new SQL({ ...this.t, schema_up: sql.join('\n'), })
  }
  down(...sql: string[]) {
    return new SQL({ ...this.t, schema_down: sql.join('\n'), })
  }
  build() {
    const { schema_up, schema_down } = this.t
    return { schema_up, schema_down }
  }
}
