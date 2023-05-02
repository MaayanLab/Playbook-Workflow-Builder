import { MaybeArray } from '@/utils/types'
export type Icon = MaybeArray<{ path: string, title?: string, transform?: string, size?: number }>
export * from './mdi'
export * from './services'