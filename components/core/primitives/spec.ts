import type { Icon } from '@/icons'
import type { InternalIdentifiableMetaNode } from '@/spec/metanode'

export type Primative = {
  name: string,
  label: string,
  icon?: Icon,
  color?: string,
  extra?: {
    term?: {
      meta?: Partial<InternalIdentifiableMetaNode['meta']>,
      autocomplete?: (search: string) => { items: string[], error?: string }
    },
    set?: {
      meta?: Partial<InternalIdentifiableMetaNode['meta']>,
    },
    scored?: {
      meta?: Partial<InternalIdentifiableMetaNode['meta']>,
    },
  }
}
