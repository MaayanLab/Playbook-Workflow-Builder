import { Primative } from './spec'
import useSWRImmutable from 'swr/immutable'
import levenSort from '@/utils/leven-sort'
import fetcher from '@/utils/next-rest-fetcher'
import { glygen_icon } from '@/icons'

// XXX
// function useHarmonizomeGeneSuggestions(search: string) {
//   const { data, error } = useSWRImmutable<string[]>(() => search.length >= 2 ? `https://maayanlab.cloud/Harmonizome/api/1.0/suggest?t=protein&q=${encodeURIComponent(search)}` : null, fetcher)
//   const items = data ? levenSort(data, search).slice(0, 10) as string[] : []
//   return { items, error }
// }

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
      // autocomplete: useHarmonizomeGeneSuggestions,
    },
    set: {
      meta: {
        example: { description: 'Example glycan set', set: ['G08146BT', 'G86750HK', 'G93822DY', 'G81459MV', 'G37626XQ',] },
        pagerank: 6,
      },
    },
  },
} as Primative
