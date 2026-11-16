import { Primative } from './spec'
import { useQuery } from '@tanstack/react-query'
import levenSort from '@/utils/leven-sort'
import { drug_icon } from '@/icons'

function usePubchemDrugSuggestions(search: string) {
  const { data, error } = useQuery({
    queryKey: ['pubchem-drug-suggestions', search],
    enabled: search.length >= 3,
    queryFn: async () => {
      const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(search)}`)
      if (!response.ok) throw new Error(await response.text())
      const { dictionary_terms: { compound } } = await response.json()
      return compound as string[]
    },
    staleTime: Infinity,
  })
  const items = data ? levenSort(data, search).slice(0, 10) as string[] : []
  return { items, error: error?.message }
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
        example: { description: 'Example drug set', set: ['ac1ndss5', 'adoprazine', 'ai-10-49', 'alisporivir', 'almitrine', 'alvocidib', 'am 580', 'amg-9810', 'amuvatinib', 'amuvatinib', 'antimycin a', 'apixaban', 'as-252424', 'avasimibe', 'avatrombopag', 'bp-897', 'brexpiprazole', 'brivanib', 'camostat', 'carboxyamidotriazole', 'cbipes', 'cc-223', 'cetylpyridinium chloride', 'chlormidazole', 'ci-1040', 'cloconazole', 'convallatoxin', 'cycloheximide', 'cyclopiazonic acid', 'cypermethrin', 'dapivirine', 'dcpib', 'deguelin', 'digoxigenin', 'dihydromunduletone', 'dihydrorotenone', 'diydroxyflavone', 'drotaverine', 'ethaverine', 'etifoxine', 'fenretinide', 'flunarizine', 'gedunin', 'gitoxigenin diacetate', 'gsk2606414', 'harringtonine', 'hematoporphyrin', 'homoharringtonine', 'imd0354', 'ipag', 'isorotenone', 'jte-013', 'ketoconazole', 'lde225', 'leoidin', 'lgk-974', 'lidoflazine', 'lonafarnib', 'lopinavir', 'loratadine', 'loteprednol etabonate', 'ly2228820', 'mefloquine', 'methylene blue', 'mibampator', 'mk-886', 'mundulone', 'nafamostat', 'nsc319726', 'octenidine', 'oxiconazole', 'papaverine', 'pevonedistat', 'pexidartinib', 'pf-670462', 'ph-797804', 'polidocanol', 'posaconazole', 'proscillaridin', 'raf265 derivative', 'ravuconazole', 'regorafenib', 'sb-612111', 'silmitasertib', 'sorafenib', 'stf-62247', 'strophanthidin', 'strophanthidinic acid', 'thapsigargin', 'thimerosal', 'thioguanosine', 'tioguanine', 'torin 1', 'torin 2', 'tyrphostin', 'vatalanib', 'vlx600', 'voxtalisib', 'vu 0155069', 'way-600', 'zk-93423'] },
      },
    },
  },
} as Primative
