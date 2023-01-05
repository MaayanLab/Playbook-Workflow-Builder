/**
 * A key-value database
 */
import levelup, { LevelUp } from 'levelup'
import memdown from 'memdown'
import sqldown from 'sqldown'

declare global {
  var kvdb: LevelUp
}
const kvdb = global.kvdb ? global.kvdb : levelup(process.env.DATABASE_URL ? sqldown(process.env.DATABASE_URL) : memdown())
global.kvdb = kvdb
export default kvdb