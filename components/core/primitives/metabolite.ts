import { Primative } from './spec'
import { drug_icon } from '@/icons'

export const Metabolite = {
  name: 'Metabolite',
  label: 'Metabolite',
  icon: [drug_icon],
  color: '#ffb997',
  extra: {
    term: {
      meta: {
        example: 'Glucose',
      },
    },
  },
} as Primative
