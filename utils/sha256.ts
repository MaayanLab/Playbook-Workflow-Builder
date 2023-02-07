import { Hash } from "fast-sha256"
import { canonicalize } from 'json-canonicalize'

/**
 * Obtain a unique hash for arbitrary data
 */
export default function sha256(data: any) {
  const h = new Hash()
  h.update(Buffer.from(canonicalize(data)))
  return Buffer.from(h.digest()).toString('hex')
}
