import { Primative } from './spec'
import { tissue_icon } from '@/icons'

export const TumorGeneExpression = {
  name: 'TumorGeneExpression',
  label: 'Tumor Gene Expression',
  icon: [tissue_icon],
  color: '#98D7C2',
  extra: {
    term: {
      meta: {
        example: 'Meningioma',
      },
    },
  },
} as Primative
