import { MetaNode } from '@/spec/metanode'
import { kea3_icon } from '@/icons'
import { z } from 'zod'
import { GeneSet } from '@/components/core/set'
import { GeneRanked } from '@/components/core/ranked'
import Citable from '@/utils/citations'

const kea3_url = 'https://maayanlab.cloud/kea3'

const KEA3ResponseT = z.record(
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

// TODO: show all KEA results, allow extracting some specific kinases
// export const KEA3KinaseEnrichmentAnalysisResults = MetaNode('KEA3KinaseEnrichmentAnalysisResults')
//   .meta({
//     label: 'Kinase Enrichment Analysis Results',
//     description: '',
//     icon: [kea3_icon, weighted_icon],
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

export const KEA3KinaseEnrichmentAnalysis = MetaNode('KEA3KinaseEnrichmentAnalysis')
  .meta({
    label: 'Kinase Enrichment Analysis',
    description: 'Use KEA3 API to infer upstream kinases whose putative substrates are overrepresented in the input set',
    icon: [kea3_icon],
  })
  .inputs({ gene_set: GeneSet })
  .output(GeneRanked)
  .resolve(async (props) => {
    const query_name = `pwb ${props.inputs.gene_set.description || 'gene set'}`
    const gene_set = props.inputs.gene_set.set
    const req = await fetch(
      `${kea3_url}/api/enrich/`,
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
    const { ['Integrated--meanRank']: results } = KEA3ResponseT.parse(await req.json())
    results.sort((a, b) => (+a.Rank) - (+b.Rank))
    return { ranked: results.map(({ TF }) => TF) }
  })
  .story(props => ({
    abstract: Citable.text`Kinase Enrichment Analysis was performed on ${props.inputs?.gene_set?.description ? props.inputs.gene_set.description : 'the gene set'} using KEA3 [${Citable.doi('10.1093/nar/gkab359')}].`
  }))
  .build()
