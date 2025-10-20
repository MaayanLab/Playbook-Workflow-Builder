import { z } from 'zod'
import { Codec } from '@/spec/codec'

function identity<T>(value: T): T {
  return value
}

export default function codecFrom<Output = any, Def extends z.ZodTypeDef = z.ZodTypeDef>(typ: z.ZodType<Output, Def, Output>): Codec<Output, Output> {
  return {
    encode: (obj: Output) => identity(obj),
    decode: (obj: Output) => typ.parse(identity(obj)),
    zod: typ,
  }
}
