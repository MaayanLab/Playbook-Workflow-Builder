/**
 * A lightweight "ORM" for constructing type-safe pg tables, useful so
 *  we can define the .sql and the typescript model at the same time.
 **/
import { z } from "zod"
import * as dict from '@/utils/dict'

export type TypedSchema<T = {}> = z.ZodType<T> & {
  name: string
  field_types: {[K in keyof T]: z.ZodType<T[K]>}
  field_sql: {[K in keyof T]: string}
  extra_insert_sql?: string,
  schema_up: string
  schema_down: string
}

export type TableSchema<T = {}> = {
  name: string
  field_types: {[K in keyof T]: z.ZodType<T[K]>}
  field_sql: {[K in keyof T]: string}
  field_extra_sql: {[K in keyof T]: string}
  extra_sql: string[]
  extra_insert_sql: string[],
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
    return new Table({ name, field_types: {}, field_sql: {}, field_extra_sql: {}, extra_sql: [], extra_insert_sql: [] })
  }
  field<S extends string, C = undefined>(name: S, sql: string, extra_sql: string, type: z.ZodType<C>) {
    return new Table({
      ...this.t,
      field_types: { ...this.t.field_types, [name]: type },
      field_sql: { ...this.t.field_sql, [name]: sql },
      field_extra_sql: { ...this.t.field_sql, [name]: extra_sql },
    } as TableSchema<T & { [K in S]: C }>)
  }
  extra(extra_sql: string) {
    return new Table({
      ...this.t,
      extra_sql: [...this.t.extra_sql, extra_sql],
    })
  }
  extra_insert(extra_insert_sql: string) {
    return new Table({
      ...this.t,
      extra_insert_sql: [...this.t.extra_insert_sql, extra_insert_sql],
    })
  }
  build() {
    const { name, field_types, extra_insert_sql } = this.t
    const type = z.object(this.t.field_types)
    const schema_up = [
      `create table ${this.t.name} (`,
      [
        ...dict.keys(this.t.field_sql).map(field =>
          `  ${JSON.stringify(field)} ${this.t.field_sql[field]} ${this.t.field_extra_sql[field]}`
        ),
        ...this.t.extra_sql.map(sql => `  ${sql}`),
      ].join(',\n'),
      `);`,
    ].join('\n')
    
    const schema_down = `drop table ${this.t.name};`
    return { ...type, name, field_types, extra_insert_sql, schema_up, schema_down } as unknown as TypedSchema<T>
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
    const { name, field_types } = this.t
    const type = z.object(this.t.field_types)
    const schema_up = `create view ${this.t.name} as ${this.t.sql}`
    const schema_down = `drop view ${this.t.name};`
    return { ...type, name, field_types, schema_up, schema_down } as unknown as TypedSchema<T>
  }
}
