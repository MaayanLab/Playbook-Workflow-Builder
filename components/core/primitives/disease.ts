import { Primative } from './spec'
import { disease_icon } from '@/icons'

export const Disease = {
  name: 'Disease',
  label: 'Disease',
  icon: [disease_icon],
  color: '#B4F8C8',
  extra: {
    term: {
      meta: {
        example: 'Diabetic Nephropathy',
      },
    },
    set: {
      meta: {
        example: { description: 'Example disease set', set: ['Diabetic Nephropathy'] },
      },
    },
  },
} as Primative
