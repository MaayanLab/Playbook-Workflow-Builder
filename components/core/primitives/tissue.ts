import { Primative } from './spec'
import { tissue_icon } from '@/icons'

export const Tissue = {
  name: 'Tissue',
  label: 'Tissue',
  icon: [tissue_icon],
  color: '#FFC8DD',
  extra: {
    term: {
      meta: {
        example: 'Brain',
      },
    },
    set: {
      meta: {
        example: { description: 'Example tissue set', set: ['Brain']},
      }
    },
  },
} as Primative
