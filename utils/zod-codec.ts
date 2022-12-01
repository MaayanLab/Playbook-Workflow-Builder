import { z } from 'zod'
import { Codec } from '@/spec/codec'

export default function codecFrom<T>(typ: z.ZodType<T>): Codec<T> {
  return {
    encode: (obj: T) => JSON.stringify(obj),
    decode: (obj: string) => typ.parse(JSON.parse(obj)),
  }
}
