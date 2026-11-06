import { z } from 'zod'

export const z_uuid = z.string
export const z_maybe_array = <Output, Def extends z.ZodTypeDef = z.ZodTypeDef, Input = Output>(inner: z.ZodType<Output, Def, Input>) => z.union([inner, z.array(inner)])
