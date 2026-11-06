import { z } from 'zod'

export const z_uuid = z.string
export const z_maybe_array = <Schema extends z.ZodType>(inner: Schema) => z.union([inner, z.array(inner)])
