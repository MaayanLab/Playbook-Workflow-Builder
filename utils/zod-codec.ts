import { z } from 'zod'
import { Codec } from '@/spec/codec'

function identity<T>(value: T): T {
  return value
}

export default function codecFrom<Output = any>(typ: z.ZodType<Output>): Codec<Output, unknown> {
  return {
    encode: (obj: Output) => identity(obj),
    decode: (obj: unknown) => typ.parse(identity(obj)),
    zod: typ,
  }
}
