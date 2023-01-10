import * as db from '@/db'
import { TypedSchema } from '@/spec/sql'

console.log([
  '-- migrate:up',
  ...Object.values<TypedSchema>(db as any).map(element => element.schema_up),
  '',
  '-- migrate:down',
  ...Object.values<TypedSchema>(db as any).reverse().map(element => element.schema_down),
].join('\n'))
