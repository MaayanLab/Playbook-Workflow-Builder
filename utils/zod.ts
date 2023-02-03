import { z } from 'zod'

export const z_uuid = z.string
export const nullable_undefined_codec = <C>(type: z.ZodType<C>) => ({
  decode: type.nullable().transform(v => v !== null ? v : undefined).parse,
  encode: type.optional().transform(v => v !== undefined ? v : null).parse,
})
