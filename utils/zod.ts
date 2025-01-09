import { z } from 'zod'

export const z_uuid = z.string
export const nullable_undefined_codec = <Output, Def extends z.ZodTypeDef = z.ZodTypeDef, Input = Output>(type: z.ZodType<Output, Def, Input>) => ({
  decode: type.nullable().transform(v => v !== null ? v : undefined).parse,
  encode: type.optional().transform(v => v !== undefined ? v : null).parse,
})
export const json_safe_timestamp_nullable_codec = () => ({
  encode: z.union([z.null(), z.date().transform((v) => v.toISOString())]).parse,
  decode: z.union([z.null(), z.date(), z.string().transform((v) => new Date(v))]).parse,
})
export const json_safe_timestamp_codec = () => ({
  encode: z.date().transform((v) => v.toISOString()).parse,
  decode: z.union([z.date(), z.string().transform((v) => new Date(v))]).parse,
})
export const z_bigint_codec = () => ({
  decode: (db: unknown) => z.union([z.number(), z.string()]).transform(v => +v).parse(db),
  encode: (v: number) => v.toString(),
})

export const z_maybe_array = <Output, Def extends z.ZodTypeDef = z.ZodTypeDef, Input = Output>(inner: z.ZodType<Output, Def, Input>) => z.union([inner, z.array(inner)])
