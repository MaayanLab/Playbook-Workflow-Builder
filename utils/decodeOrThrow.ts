import * as t from 'io-ts'
import * as Either from 'fp-ts/Either'
import { PathReporter } from 'io-ts/lib/PathReporter'

export default function decodeOrThrow<A, O, I>(typ: t.Type<A, O, I>, i: I): A {
  return Either.getOrElse<t.Errors, A>(
    (err) => { throw new Error(PathReporter.report(Either.left(err)).join('\n')) }
  )(typ.decode(i))
}
