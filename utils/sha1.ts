import crypto from 'crypto'

/**
 * Obtain a unique hash for a arbitrary data
 */
export default function sha1(data: any) {
  const h = crypto.createHash('sha1')
  h.update(Buffer.from(JSON.stringify(data)))
  return h.digest('hex')
}

export function sha1_to_uuid(sha1: string) {
  return [
    sha1.substring(0, 8),
    sha1.substring(8, 12),
    sha1.substring(12, 16),
    sha1.substring(16, 20),
    sha1.substring(20, 32),
  ].join('-')
}

export function uuid_to_sha1(uuid: string) {
  return uuid.replaceAll(/-/g, '')
}
