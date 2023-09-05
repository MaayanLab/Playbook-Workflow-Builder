import { MetaNode } from '@/spec/metanode'
import { ScoredTissues } from '@/components/core/input/scored'
import { archs4_icon } from '@/icons'
import { GeneTerm } from '@/components/core/input/term'
import read_csv from '@/utils/csv'

async function archs4_tissue_expression({ search, species = 'human', type = 'tissue' }: { search: string, species?: string, type?: string }) {
  const params = new URLSearchParams()
  params.append('search', search)
  params.append('species', species)
  params.append('type', type)
  const req = await fetch(`https://maayanlab.cloud/archs4/search/loadExpressionTissue.php?${params.toString()}`)
  const { values } = read_csv<'id' | 'min' | 'q1' | 'median' | 'q3' | 'max' | 'color'>(await req.text())
  const results = values
    .map(({ id, median }) => ({
      term: id.split('.')[3],
      zscore: (Number(median) - 4.336)/3.921, // these mu/std coefficients were obtained by querying random genes
    }))
    .filter(({ term }) => !!term)
  results.sort((a, b) => b.zscore - a.zscore)
  return results
}

export const ARCHS4TissueExpression = MetaNode('ARCHS4TissueExpression')
  .meta({
    label: 'Query ARCHS4 Median Tissue Expression',
    description: 'Use ARCHS4 API to obtain median tissue expression for the given gene',
    icon: [archs4_icon],
  })
  .inputs({ gene: GeneTerm })
  .output(ScoredTissues)
  .resolve(async (props) => {
    return await archs4_tissue_expression({ search: props.inputs.gene })
  })
  .story(props =>
    `Median expression of ${props.inputs ? props.inputs.gene : 'the gene'} was obtained from ARCHS4 [\\ref{doi:10.1038/s41467-018-03751-6}].`
  )
  .build()
