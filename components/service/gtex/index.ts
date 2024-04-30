import { MetaNode } from '@/spec/metanode'
import { GeneInfo, GeneInfoFromGeneTerm } from '@/components/service/mygeneinfo'
import { ScoredTissues } from '@/components/core/scored'
import { gtex_icon } from '@/icons'
import python from '@/utils/python'
import { GeneTerm } from '@/components/core/term'

export const GTExTissueExpression = MetaNode('GTExTissueExpression')
  .meta({
    label: 'Query GTEx Median Tissue Expression',
    description: 'Use GTEx API to obtain median tissue expression for the given gene',
    icon: [gtex_icon],
    pagerank: 1,
  })
  .inputs({ gene_info: GeneInfo })
  .output(ScoredTissues)
  .resolve(async (props) => {
    return await python(
      'components.service.gtex.gtex_gene_expression',
      { kargs: [props.inputs.gene_info.ensembl?.gene || props.inputs.gene_info.symbol], kwargs: { datasetId: 'gtex_v8' } },
      message => props.notify({ type: 'info', message }),
    )
  })
  .story(props =>
    `Median expression of ${props.inputs ? props.inputs.gene_info.symbol : 'the gene'} was obtained from the GTEx Portal [\\ref{doi:10.1038/ng.2653}] using the portal's API.`
  )
  .build()

export const GTExTissueExpressionFromGene = MetaNode('GTExTissueExpressionFromGene')
  .meta(GTExTissueExpression.meta)
  .inputs({ gene: GeneTerm })
  .output(GTExTissueExpression.output)
  .resolve(async (props) => {
    const gene_info = await GeneInfoFromGeneTerm.resolve(props)
    return await GTExTissueExpression.resolve({ ...props, inputs: { gene_info } })
  })
  .story(props =>
    `Median expression of ${props.inputs ? props.inputs.gene : 'the gene'} was obtained from the GTEx Portal [\\ref{doi:10.1038/ng.2653}] using the portal's API.`
  )
  .build()
