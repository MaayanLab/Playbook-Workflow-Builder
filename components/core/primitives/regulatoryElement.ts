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
        example: { description: 'Example regulatory-element set', set: ['EH38E2744005','EH38E2744007','EH38E2744008','EH38E3926310','EH38E3926311','EH38E2744014','EH38E3926315','EH38E3926316']},
        pagerank: 6,
      }
    },
  },
} as Primative
