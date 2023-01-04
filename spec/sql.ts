/**
 * A lightweight "ORM" for constructing type-safe pg tables, useful so
 *  we can define the .sql and the typescript model at the same time.
 **/
import { z } from "zod"
import * as dict from '@/utils/dict'

export type TypedSchema<T = {}> = z.ZodType<T> & {
  schema_up: string
  schema_down: string
}

export type TableSchema<T = {}> = {
  name: string
  field_types: {[K in keyof T]: z.ZodType<T[K]>}
  field_sql: {[K in keyof T]: string}
  extra_sql: string[]
}

/**
 * Usage:
 * const table_name = Table.createTable('table_name')
 *   .field('fieldname', z.string(), 'varchar') // fieldname, zod schema and pg schema
 *   .extra('') // any extra sql in the table definition like primary key/constraints
 *   .build()
 */
export class Table<T = {}> {
  constructor(public t: TableSchema<T>) {}
  static create(name: string) {
    return new Table({ name, field_types: {}, field_sql: {}, extra_sql: [] })
  }
  field<S extends string, C = undefined>(name: S, sql: string, type: z.ZodType<C>) {
    return new Table({
      ...this.t,
      field_types: { ...this.t.field_types, [name]: type },
      field_sql: { ...this.t.field_sql, [name]: sql },
    } as TableSchema<T & { [K in S]: C }>)
  }
  extra(sql: string) {
    return new Table({
      ...this.t,
      extra_sql: [...this.t.extra_sql, sql],
    })
  }
  build() {
    const type = z.object(this.t.field_types)
    const schema_up = [
      `create table ${this.t.name} (`,
      [
        ...dict.items(this.t.field_sql).map(
          ({ key, value: sql }) => `  ${JSON.stringify(key)} ${sql}`
        ),
        ...this.t.extra_sql.map(sql => `  ${sql}`),
      ].join(',\n'),
      `);`,
    ].join('\n')
    
    const schema_down = `drop table ${this.t.name};`
    return { ...type, schema_up, schema_down } as unknown as TypedSchema<T>
  }
}

export type ViewSchema<T = {}> = {
  name: string
  field_types: {[K in keyof T]: z.ZodType<T[K]>}
  sql: string
}

export class View<T extends {} = {}> {
  constructor(public t: ViewSchema<T>) {}
  static create(name: string) {
    return new View({ name, field_types: {}, sql: '' })
  }
  field<S extends string, C = undefined>(name: S, type: z.ZodType<C>) {
    return new View({
      ...this.t,
      field_types: { ...this.t.field_types, [name]: type },
    } as ViewSchema<T & { [K in S]: C }>)
  }
  sql(...sql: string[]) {
    return new View({ ...this.t, sql: sql.join('\n'), })
  }
  build() {
    const type = z.object(this.t.field_types)
    const schema_up = `create view ${this.t.name} as ${this.t.sql}`
    const schema_down = `drop view ${this.t.name};`
    return { ...type, schema_up, schema_down } as unknown as TypedSchema<T>
  }
}
