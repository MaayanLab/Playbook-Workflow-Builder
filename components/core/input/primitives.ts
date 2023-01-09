import useSWRImmutable from 'swr/immutable'
import levenSort from '@/utils/leven-sort'
import { gene_icon, drug_icon } from '@/icons'

const fetcher = (url: string) => fetch(url).then(r => r.json())

function useHarmonizomeGeneSuggestions(search: string) {
  const { data, error } = useSWRImmutable<string[]>(() => search.length >= 2 ? `https://maayanlab.cloud/Harmonizome/api/1.0/suggest?t=gene&q=${encodeURIComponent(search)}` : null, fetcher)
  const items = data ? levenSort(data, search).slice(0, 10) as string[] : []
  return { items, error }
}

export const Gene = {
  name: 'Gene',
  label: 'Gene',
  icon: [gene_icon],
  color: '#B3CFFF',
  examples: {
    term: 'ACE2',
    set: ['UTP14A', 'S100A6', 'SCAND1', 'RRP12', 'CIAPIN1', 'ADH5', 'MTERF3', 'SPR', 'CHMP4A', 'UFM1', 'VAT1', 'HACD3', 'RFC5', 'COTL1', 'NPRL2', 'TRIB3', 'PCCB', 'TLE1', 'CD58', 'BACE2', 'KDM3A', 'TARBP1', 'RNH1', 'CHAC1', 'MBNL2', 'VDAC1', 'TES', 'OXA1L', 'NOP56', 'HAT1', 'CPNE3', 'DNMT1', 'ARHGAP1', 'VPS28', 'EIF2S2', 'BAG3', 'CDCA4', 'NPDC1', 'RPS6KA1', 'FIS1', 'SYPL1', 'SARS', 'CDC45', 'CANT1', 'HERPUD1', 'SORBS3', 'MRPS2', 'TOR1A', 'TNIP1', 'SLC25A46', 'MAL', 'EPCAM', 'HDAC6', 'CAPN1', 'TNRC6B', 'PKD1', 'RRS1', 'HP', 'ANO10', 'CEP170B', 'IDE', 'DENND2D', 'CAMK2B', 'ZNF358', 'RPP38', 'MRPL19', 'NUCB2', 'GNAI1', 'LSR', 'ADGRE2', 'PKMYT1', 'CDK5R1', 'ABL1', 'PILRB', 'AXIN1', 'FBXL8', 'MCF2L', 'DBNDD1', 'IGHMBP2', 'WIPF2', 'WFS1', 'OGFOD2', 'MAPK1IP1L', 'COL11A1', 'REG3A', 'SERPINA1', 'MYCBP2', 'PIGK', 'TCAP', 'CRADD', 'ELK1', 'DNAJB2', 'ZBTB16', 'DAZAP1', 'MAPKAPK2', 'EDRF1', 'CRIP1', 'UCP3', 'AGR2', 'P4HA2',],
  },
  autocomplete: {
    term: useHarmonizomeGeneSuggestions,
  },
}

const pubchemFetcher = (url: string) => fetcher(url).then(({ dictionary_terms: { compound } }) => compound)

function usePubchemDrugSuggestions(search: string) {
  const { data, error } = useSWRImmutable<string[]>(() => search.length >= 3 ? `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(search)}` : null, pubchemFetcher)
  const items = data ? levenSort(data, search).slice(0, 10) as string[] : []
  return { items, error }
}

export const Drug = {
  name: 'Drug',
  label: 'Drug',
  icon: [drug_icon],
  color: '#FFE4A0',
  examples: {
    term: 'imatinib',
    set: ['ac1ndss5', 'adoprazine', 'ai-10-49', 'alisporivir', 'almitrine', 'alvocidib', 'am 580', 'amg-9810', 'amuvatinib', 'amuvatinib', 'antimycin a', 'apixaban', 'as-252424', 'avasimibe', 'avatrombopag', 'bp-897', 'brexpiprazole', 'brivanib', 'camostat', 'carboxyamidotriazole', 'cbipes', 'cc-223', 'cetylpyridinium chloride', 'chlormidazole', 'ci-1040', 'cloconazole', 'convallatoxin', 'cycloheximide', 'cyclopiazonic acid', 'cypermethrin', 'dapivirine', 'dcpib', 'deguelin', 'digoxigenin', 'dihydromunduletone', 'dihydrorotenone', 'diydroxyflavone', 'drotaverine', 'ethaverine', 'etifoxine', 'fenretinide', 'flunarizine', 'gedunin', 'gitoxigenin diacetate', 'gsk2606414', 'harringtonine', 'hematoporphyrin', 'homoharringtonine', 'imd0354', 'ipag', 'isorotenone', 'jte-013', 'ketoconazole', 'lde225', 'leoidin', 'lgk-974', 'lidoflazine', 'lonafarnib', 'lopinavir', 'loratadine', 'loteprednol etabonate', 'ly2228820', 'mefloquine', 'methylene blue', 'mibampator', 'mk-886', 'mundulone', 'nafamostat', 'nsc319726', 'octenidine', 'oxiconazole', 'papaverine', 'pevonedistat', 'pexidartinib', 'pf-670462', 'ph-797804', 'polidocanol', 'posaconazole', 'proscillaridin', 'raf265 derivative', 'ravuconazole', 'regorafenib', 'sb-612111', 'silmitasertib', 'sorafenib', 'stf-62247', 'strophanthidin', 'strophanthidinic acid', 'thapsigargin', 'thimerosal', 'thioguanosine', 'tioguanine', 'torin 1', 'torin 2', 'tyrphostin', 'vatalanib', 'vlx600', 'voxtalisib', 'vu 0155069', 'way-600', 'zk-93423'],
  },
  autocomplete: {
    term: usePubchemDrugSuggestions,
  },
}

export type Primative = typeof Gene | typeof Drug
