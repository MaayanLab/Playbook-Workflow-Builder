import { z } from 'zod'

export type Codec<D, E = string> = {
  encode: (data: D) => E
  decode: (bytes: E) => D
  zod?: z.ZodType<D>,
}
export type Decoded<T> = T extends Codec<infer D, infer _E> ? D : never
export type Encoded<T> = T extends Codec<infer _D, infer E> ? E : never