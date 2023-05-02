import * as db from '@/db'
import { TypedSchema } from '@/spec/sql'
import * as dict from '@/utils/dict'

console.log([
  '-- migrate:up',
  ...dict.values(db as any)
    .filter((item): item is Pick<TypedSchema<unknown>, 'schema_up'> => 'schema_up' in item)
    .map((element) => element.schema_up),
  '',
  '-- migrate:down',
  ...dict.values(db as any)
    .filter((item): item is Pick<TypedSchema<unknown>, 'schema_down'> => 'schema_down' in item)
    .reverse()
    .map((element) => element.schema_down),
].join('\n'))
