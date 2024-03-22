import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/input/term'
import { z } from 'zod'
import { gene_icon, mygeneinfo_icon } from '@/icons'

export const MyGeneInfoHitC = z.object({
  hits: z.array(
    z.object({
      _id: z.string(),
      _score: z.number(),
      symbol: z.string().optional(),
      name: z.string().optional(),
      taxid: z.number().optional(),
      entrezgene: z.string().optional(),
    }),
  ),
})
export type MyGeneInfoHit = z.infer<typeof MyGeneInfoHitC>

export const MyGeneInfoC = z.object({
  _id: z.string(),
  symbol: z.string(),
  entrezgene: z.string().optional(),
  name: z.string().optional(),
  taxid: z.number().optional(),
  ensembl: z.object({
    gene: z.string(),
  }).optional()
})
export type MyGeneInfo = z.infer<typeof MyGeneInfoC>

async function mygeneinfo_query(geneId: string): Promise<{total: number, hits: Array<MyGeneInfoHit>}> {
  const res = await fetch(`https://mygene.info/v3/query?q=${encodeURIComponent(geneId)}`)
  return await res.json()
}

async function mygeneinfo(geneId: string): Promise<MyGeneInfo> {
  const res = await fetch(`https://mygene.info/v3/gene/${encodeURIComponent(geneId)}`)
  return await res.json()
}

async function getGeneData(geneSymbol: string){
  const results = MyGeneInfoHitC.parse(await mygeneinfo_query(geneSymbol))
  const hits = results.hits.filter((hit): hit is MyGeneInfoHit['hits'][0] & { symbol: string } => !!hit.symbol)
  const exactMatch = hits.filter(hit => hit.symbol.toUpperCase() == geneSymbol.toUpperCase())[0]
  const _id: string | undefined = exactMatch !== undefined ? exactMatch._id : hits[0]._id
  if (_id === undefined) throw new Error(`Could not identify a gene for the symbol ${geneSymbol} in mygene.info`)
  return await mygeneinfo(_id)
}

export const GeneInfo = MetaNode('GeneInfo')
  .meta({
    label: 'Gene Information',
    description: 'A Gene resolved with MyGeneInfo',
    icon: [gene_icon],
  })
  .codec(MyGeneInfoC)
  .view(geneinfo => (
    <div className="prose max-w-none">
      <a href={`https://www.ncbi.nlm.nih.gov/gene/${geneinfo.entrezgene}`}>{geneinfo.symbol}</a> {geneinfo.name}
    </div>
  ))
  .build()

export const GeneInfoFromGeneTerm = MetaNode('GeneInfoFromGeneTerm')
  .meta({
    label: 'Resolve Gene Info from Term',
    description: 'Resolve gene info from gene term with MyGeneInfo',
    icon: [mygeneinfo_icon],
  })
  .inputs({ gene: GeneTerm })
  .output(GeneInfo)
  .resolve(async (props) => {
    return await getGeneData(props.inputs.gene);
  })
  .story(props =>
    `More information about the gene was then obtained with the MyGene.info API [\\ref{doi:10.1186/s13059-016-0953-9},\\ref{doi:10.1093/nar/gks1114}].`
  )
  .build()
