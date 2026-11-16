import { Primative } from './spec'
import { glygenprotein_icon } from '@/icons'

export const Protein = {
  name: 'Protein',
  label: 'Protein',
  icon: [glygenprotein_icon],
  color: '#B3CFFF',
  extra: {
    term: {
      meta: {
        example: 'P14210',
        pagerank: 7,
      },
    },
    set: {
      meta: {
        example: { description: 'Example protein set', set: ['P14210', 'P29274',] },
        pagerank: 6,
      },
    },
  },
} as Primative
