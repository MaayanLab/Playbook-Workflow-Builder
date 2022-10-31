/**
 * A key-value database
 */
import levelup, { LevelUp } from 'levelup'
import leveldown from 'leveldown'
declare global {
  var kvdb: LevelUp
}
const kvdb = global.kvdb ? global.kvdb : levelup(leveldown('data'))
global.kvdb = kvdb
export default kvdb