export type Codec<D, E = string> = {
  encode: (data: D) => E
  decode: (bytes: E) => D
}