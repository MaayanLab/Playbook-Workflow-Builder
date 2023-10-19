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
        example: { description: 'Example variantset', set: ['CA13203640','CA659845332','CA659845333','CA932539093','CA2525517560','CA595730759','CA213257739','CA213257743','CA213257744','CA213257745','CA213257748','CA213257754']},
        pagerank: 6,
      }
    },
  },
} as Primative
