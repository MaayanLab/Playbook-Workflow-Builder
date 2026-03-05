import { Primative } from './spec'
import { variable_icon } from '@/icons'

export const PMCAccession = {
  name: 'PMC Accession',
  label: 'PMC Accession',
  icon: [variable_icon],
  color: '#8dd3c7',
  extra: {
    term: {
      meta: {
        example: 'PMC12881680',
      },
    },
  },
} as Primative
