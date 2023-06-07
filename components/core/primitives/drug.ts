import { Primative } from './spec'
import useSWRImmutable from 'swr/immutable'
import levenSort from '@/utils/leven-sort'
import fetcher from '@/utils/next-rest-fetcher'
import { drug_icon } from '@/icons'

const pubchemFetcher = (url: string) => fetcher<any>(url).then(({ dictionary_terms: { compound } }) => compound)

function usePubchemDrugSuggestions(search: string) {
  const { data, error } = useSWRImmutable<string[]>(() => search.length >= 3 ? `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(search)}` : null, pubchemFetcher)
  const items = data ? levenSort(data, search).slice(0, 10) as string[] : []
  return { items, error }
}

export const Drug = {
  name: 'Drug',
  label: 'Drug',
  icon: [drug_icon],
  color: '#FBE7C6',
  extra: {
    term: {
      meta: {
        example: 'imatinib',
        pagerank: 3,
      },
      autocomplete: usePubchemDrugSuggestions,
    },
    set: {
      meta: {
        example: { description: 'Example drugset', set: ['ac1ndss5', 'adoprazine', 'ai-10-49', 'alisporivir', 'almitrine', 'alvocidib', 'am 580', 'amg-9810', 'amuvatinib', 'amuvatinib', 'antimycin a', 'apixaban', 'as-252424', 'avasimibe', 'avatrombopag', 'bp-897', 'brexpiprazole', 'brivanib', 'camostat', 'carboxyamidotriazole', 'cbipes', 'cc-223', 'cetylpyridinium chloride', 'chlormidazole', 'ci-1040', 'cloconazole', 'convallatoxin', 'cycloheximide', 'cyclopiazonic acid', 'cypermethrin', 'dapivirine', 'dcpib', 'deguelin', 'digoxigenin', 'dihydromunduletone', 'dihydrorotenone', 'diydroxyflavone', 'drotaverine', 'ethaverine', 'etifoxine', 'fenretinide', 'flunarizine', 'gedunin', 'gitoxigenin diacetate', 'gsk2606414', 'harringtonine', 'hematoporphyrin', 'homoharringtonine', 'imd0354', 'ipag', 'isorotenone', 'jte-013', 'ketoconazole', 'lde225', 'leoidin', 'lgk-974', 'lidoflazine', 'lonafarnib', 'lopinavir', 'loratadine', 'loteprednol etabonate', 'ly2228820', 'mefloquine', 'methylene blue', 'mibampator', 'mk-886', 'mundulone', 'nafamostat', 'nsc319726', 'octenidine', 'oxiconazole', 'papaverine', 'pevonedistat', 'pexidartinib', 'pf-670462', 'ph-797804', 'polidocanol', 'posaconazole', 'proscillaridin', 'raf265 derivative', 'ravuconazole', 'regorafenib', 'sb-612111', 'silmitasertib', 'sorafenib', 'stf-62247', 'strophanthidin', 'strophanthidinic acid', 'thapsigargin', 'thimerosal', 'thioguanosine', 'tioguanine', 'torin 1', 'torin 2', 'tyrphostin', 'vatalanib', 'vlx600', 'voxtalisib', 'vu 0155069', 'way-600', 'zk-93423'] },
      },
    },
  },
} as Primative
