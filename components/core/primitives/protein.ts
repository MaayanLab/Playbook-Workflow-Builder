import { Primative } from './spec'
import useSWRImmutable from 'swr/immutable'
import levenSort from '@/utils/leven-sort'
import fetcher from '@/utils/next-rest-fetcher'
import { protein_icon } from '@/icons'

function useHarmonizomeGeneSuggestions(search: string) {
  const { data, error } = useSWRImmutable<string[]>(() => search.length >= 2 ? `https://maayanlab.cloud/Harmonizome/api/1.0/suggest?t=protein&q=${encodeURIComponent(search)}` : null, fetcher)
  const items = data ? levenSort(data, search).slice(0, 10) as string[] : []
  return { items, error }
}

export const Protein = {
  name: 'Protein',
  label: 'Protein',
  icon: [protein_icon],
  color: '#B3CFFF',
  extra: {
    term: {
      meta: {
        example: 'ADORA2A',
        pagerank: 7,
      },
      // autocomplete: useHarmonizomeGeneSuggestions,
    },
    set: {
      meta: {
        example: { description: 'Example protein set', set: ['ADORA1', 'ADORA2A', 'ADORA2B', 'ADORA3',] },
        pagerank: 6,
      },
    },
  },
} as Primative
