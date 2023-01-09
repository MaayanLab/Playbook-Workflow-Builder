import { Hash } from "fast-sha256"

/**
 * Obtain a unique hash for a arbitrary data
 */
export default function sha256(data: any) {
  const h = new Hash()
  h.update(Buffer.from(JSON.stringify(data)))
  return Buffer.from(h.digest()).toString('hex')
}
