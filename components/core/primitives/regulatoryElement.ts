import { Primative } from './spec'
import { gene_icon } from '@/icons'

export const RegulatoryElement = {
  name: 'Regulatory_Element',
  label: 'Regulatory Element',
  icon: [gene_icon],
  color: '#BDE0FE',
  extra: {
    term: {
      meta: {
        example: 'EH38E2924876',
      },
    },
    set: {
      meta: {
        example: { description: 'Example regulatory-element set', set: ['EH38E2924876', 'EH12E1234567']},
        pagerank: 6,
      }
    },
  },
} as Primative
