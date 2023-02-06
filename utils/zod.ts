import { z } from 'zod'

export const z_uuid = z.string
export const nullable_undefined_codec = <C>(type: z.ZodType<C>) => ({
  decode: type.nullable().transform(v => v !== null ? v : undefined).parse,
  encode: type.optional().transform(v => v !== undefined ? v : null).parse,
})
export const json_safe_timestamp_codec = () => ({
  decode: z.date().transform((v) => v.toISOString()).parse,
  encode: z.string().transform((v) => new Date(v)).parse,
})
