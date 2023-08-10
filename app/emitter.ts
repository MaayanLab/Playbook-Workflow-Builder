import cache from '@/utils/global_cache'
import { EventEmitter } from 'stream'

export default cache('emitter', () => new EventEmitter())
