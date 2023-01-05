export type Codec<D, E = string> = {
  encode: (data: D) => E
  decode: (bytes: E | unknown) => D
}
export type Decoded<T> = T extends Codec<infer D, infer _E> ? D : never
export type Encoded<T> = T extends Codec<infer _D, infer E> ? E : never