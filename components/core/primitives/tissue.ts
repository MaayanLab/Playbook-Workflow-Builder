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
  },
} as Primative
