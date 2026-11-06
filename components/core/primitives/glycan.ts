import { Primative } from './spec'
import { glygen_icon } from '@/icons'

export const Glycan = {
  name: 'glycan',
  label: 'Glycan',
  icon: [glygen_icon],
  color: '#B3CFFF',
  extra: {
    term: {
      meta: {
        example: 'G17689DH',
        pagerank: 7,
      },
    },
    set: {
      meta: {
        example: { description: 'Example glycan set', set: ['G08146BT', 'G86750HK', 'G93822DY', 'G81459MV', 'G37626XQ',] },
        pagerank: 6,
      },
    },
  },
} as Primative
