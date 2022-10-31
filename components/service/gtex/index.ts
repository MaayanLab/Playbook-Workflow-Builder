import { MetaNode } from '@/spec/metanode'
import { GeneInfo } from '@/components/service/mygeneinfo'
import { SignificantTissues } from '@/components/core/significant_tissues'
import { gtex_icon } from '@/icons'
import python from '@/utils/python'

export const GTExTissueExpression = MetaNode.createProcess('GTExTissueExpression')
  .meta({
    label: 'Query GTEx Median Tissue Expression',
    description: 'Use GTEx API to obtain median tissue expression for the given gene',
    icon: [gtex_icon],
  })
  .inputs({ gene_info: GeneInfo })
  .output(SignificantTissues)
  .resolve(async (props) => {
    return await python(
      'components.service.gtex.gtex_gene_expression',
      { kargs: [props.inputs.gene_info.ensembl?.gene || props.inputs.gene_info.symbol], kwargs: { datasetId: 'gtex_v8' } },
    )
  })
  .build()
