import { MetaNode } from '@/spec/metanode'
import { GeneInfo } from '@/components/mygeneinfo'
import { SignificantTissues } from '@/components/significant_tissues'
import python from '@/utils/python'

export const GTExTissueExpression = MetaNode.createProcess('GTExTissueExpression')
  .meta({
    label: 'Query GTEx Tissue eQTLs',
    description: 'Use GTEx API to identify tissue eQTLs associated with a given gene',
  })
  .inputs({ gene_info: GeneInfo })
  .output(SignificantTissues)
  .resolve(async (props) => {
    return await python(
      'components.gtex.gtex_gene_expression',
      { kargs: [props.inputs.gene_info.ensembl?.gene || props.inputs.gene_info.symbol], kwargs: { datasetId: 'gtex_v8' } },
    )
  })
  .build()
