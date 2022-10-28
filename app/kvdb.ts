/**
 * An in memory key-value database
 */
declare global {
  var kvdb: Record<string, string> | undefined
}
const kvdb: Record<string, string> = global.kvdb || {}
export default kvdb