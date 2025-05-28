import { MetaNode } from '@/spec/metanode'
import { chea3_icon } from '@/icons'
import { z } from 'zod'
import { GeneSet } from '@/components/core/set'
import { GeneRanked } from '@/components/core/ranked'

const chea3_url = 'https://maayanlab.cloud/chea3'

const ChEA3ResponseT = z.record(
  z.string(), z.array(z.object({
    'FDR': z.any(),
    'FET p-value': z.any(),
    'Intersect': z.any(),
    'Odds Ratio': z.any(),
    'Overlapping_Genes': z.any(),
    'Rank': z.any(),
    'Scaled Rank': z.any(),
    'Set_name': z.any(),
    'TF': z.any(),
  }))
)

// TODO: show all ChEA results, allow extracting some specific TFs
// export const ChEA3TFEnrichmentAnalysisResults = MetaNode('ChEA3TFEnrichmentAnalysisResults')
//   .meta({
//     label: 'Transcription Factor Enrichment Analysis Results',
//     description: '',
//     icon: [chea3_icon, weighted_icon],
//   })
//   .codec(z.record(
//     z.string(), z.array(z.object({
//       'FDR': z.any(),
//       'FET p-value': z.any(),
//       'Intersect': z.any(),
//       'Odds Ratio': z.any(),
//       'Overlapping_Genes': z.any(),
//       'Rank': z.any(),
//       'Scaled Rank': z.any(),
//       'Set_name': z.any(),
//       'TF': z.any(),
//     }))
//   ))
//   .view(props => {
//     return (
//       <>
//         {dict.items(props).map(({key, value: results}) =>
//           <Table
//             height={500}
//             cellRendererDependencies={[results]}
//             numRows={results.length}
//             enableGhostCells
//             enableFocusedCell
//             downloads={{
//               JSON: () => downloadBlob(new Blob([JSON.stringify(results)], { type: 'application/json;charset=utf-8' }), 'data.json'),
//               CSV: () => downloadBlob(new Blob([
//                 [
//                   'FDR,FET p-value,Intersect,Odds Ratio,Overlapping_Genes,Rank,Scaled Rank,Set_name,TF',
//                   ...results.map(result => `${result['FDR']},${result['FET p-value']},${result['Intersect']},${result['Odds Ratio']},${result['Overlapping_Genes']},${result['Rank']},${result['Scaled Rank']},${result['Set_name']},${result['TF']}`)
//                 ].join('\n')
//               ], { type: 'text/csv;charset=utf-8' }), 'data.csv'),
//             }}
//           >
//             <Column
//               name="FDR"
//               cellRenderer={row => <Cell key={row+''}>{results[row]['FDR']}</Cell>} />
//             <Column
//               name="FET p-value"
//               cellRenderer={row => <Cell key={row+''}>{results[row]['FET p-value']}</Cell>} />
//             <Column
//               name="Intersect"
//               cellRenderer={row => <Cell key={row+''}>{results[row]['Intersect']}</Cell>} />
//             <Column
//               name="Odds Ratio"
//               cellRenderer={row => <Cell key={row+''}>{results[row]['Odds Ratio']}</Cell>} />
//             <Column
//               name="Overlapping_Genes"
//               cellRenderer={row => <Cell key={row+''}>{results[row]['Overlapping_Genes']}</Cell>} />
//             <Column
//               name="Rank"
//               cellRenderer={row => <Cell key={row+''}>{results[row]['Rank']}</Cell>} />
//             <Column
//               name="Scaled Rank"
//               cellRenderer={row => <Cell key={row+''}>{results[row]['Scaled Rank']}</Cell>} />
//             <Column
//               name="Set_name"
//               cellRenderer={row => <Cell key={row+''}>{results[row]['Set_name']}</Cell>} />
//             <Column
//               name="TF"
//               cellRenderer={row => <Cell key={row+''}>{results[row]['TF']}</Cell>} />
//           </Table>
//         )}
//       </>
//     )
//   })
//   .build()

export const ChEA3TFEnrichmentAnalysis = MetaNode('ChEA3TFEnrichmentAnalysis')
  .meta({
    label: 'Transcription Factor Enrichment Analysis',
    description: 'Use ChEA3 API to identify TFs associated with the input set',
    icon: [chea3_icon],
    external: true,
  })
  .inputs({ gene_set: GeneSet })
  .output(GeneRanked)
  .resolve(async (props) => {
    const query_name = `pwb ${props.inputs.gene_set.description || 'gene set'}`
    const gene_set = props.inputs.gene_set.set
    const req = await fetch(
      `${chea3_url}/api/enrich/`,
      {
        method: 'POST',
        headers: {
          // 'Content-Type': 'application/json',
        },
        body: JSON.stringify({ query_name, gene_set }),
      }
    )
    if (req.status !== 200) {
      throw new Error(`Request failed (${req.status}): ${await req.text()}`)
    }
    const { ['Integrated--meanRank']: results } = ChEA3ResponseT.parse(await req.json())
    results.sort((a, b) => (+a.Rank) - (+b.Rank))
    return { ranked: results.map(({ TF }) => TF) }
  })
  .story(props => ({
    abstract: `Transcription Factor Enrichment Analysis was performed on ${props.inputs?.gene_set ? props.inputs.gene_set.description : 'the gene set'} using ChEA3\\ref{doi:10.1093/nar/gkz446}.`,
    introduction: `Transcription factors (TFs) play an important role in gene regulatory networks. TF Enrichment Analysis can be used to identify relevant TFs based on known TF-target interactions.`,
    methods: `Transcription Factor Enrichment Analysis was performed on the gene set using the ChEA3\\ref{doi:10.1093/nar/gkz446} web API and using Mean Rank across all available TF Libraries.`,
    tableLegend: `A table of significantly enriched TFs identified with ChEA3\\ref{doi:10.1093/nar/gkz446}.`,
  }))
  .build()
