import { Primative } from './spec'
import { gene_icon } from '@/icons'

export const Variant = {
  name: 'Variant',
  label: 'Variant',
  icon: [gene_icon],
  color: '#a1e8cc',
  extra: {
    term: {
      meta: {
        example: 'CA13203640',
      },
    },
    set: {
      meta: {
        example: { description: 'Example variantset', set: ['CA13203640', 'CA321211']},
        pagerank: 6,
      }
    },
  },
} as Primative
