import * as t from 'io-ts'
import { Codec } from '@/spec/codec'
import decodeOrThrow from '@/utils/decodeOrThrow'

export default function codecFrom<A, O, I>(typ: t.Type<A, O, I>): Codec<A> {
  return {
    encode: (obj: A) => JSON.stringify(obj),
    decode: (obj: string) => decodeOrThrow(typ, JSON.parse(obj)),
  }
}
