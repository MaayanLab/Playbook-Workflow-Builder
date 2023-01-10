import crypto from 'crypto'
import { canonicalize } from 'json-canonicalize'

/**
 * Obtain a unique hash for a arbitrary data
 */
export default function uuid5(data: any) {
  const h = crypto.createHash('sha1')
  h.update(Buffer.from(canonicalize(data)))
  const hex = h.digest('hex').substring(0, 32)
  return [
    hex.substring(0, 8),
    hex.substring(8, 12),
    hex.substring(12, 16),
    hex.substring(16, 20),
    hex.substring(20, 32),
  ].join('-')
}
