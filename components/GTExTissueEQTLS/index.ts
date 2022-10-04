import { MetaNode } from '@/spec/metanode'
import { GeneInfo } from '@/components/GeneInfo'
import { SignificantTissues } from '@/components/SignificantTissues'
import python from '@/utils/python'

export const GTExTissueEQTLS = MetaNode.createProcess('GTExTissueEQTLS')
  .meta({
    label: 'Query GTEx Tissue eQTLs',
    description: 'Use GTEx API to identify tissue eQTLs associated with a given gene',
  })
  .inputs({ gene_info: GeneInfo })
  .output(SignificantTissues)
  .resolve(async (props) => {
    return await python(
      '@/components/GTExTissueEQTLS/gtex_tissue.py',
      'gtex_singleTissueEqtl',
      { kargs: [props.inputs.gene_info.symbol], kwargs: { datasetId: 'gtex_v8' } },
    )
  })
  .build()
