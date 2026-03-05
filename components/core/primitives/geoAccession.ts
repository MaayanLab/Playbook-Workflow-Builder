import { Primative } from './spec'
import { variable_icon } from '@/icons'

export const GEOAccession = {
  name: 'GEO Accession',
  label: 'GEO Accession',
  icon: [variable_icon],
  color: '#8dd3c7',
  extra: {
    term: {
      meta: {
        example: 'GSE301503',
      },
    },
  },
} as Primative
