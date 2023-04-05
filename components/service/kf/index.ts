import { MetaNode } from '@/spec/metanode'
import { GeneInfo } from '@/components/service/mygeneinfo'
import { GeneExpressionInTumor } from '@/components/core/input/tumor'
import python from '@/utils/python'

export const KFTumorExpression = MetaNode('KFTumorExpression')
  .meta({
    label: 'Query KF Gene Expression in Tumor',
    description: 'Use KF API to obtain tumors expressing the given gene'
  })
  .inputs({ gene_info: GeneInfo })
  .output(GeneExpressionInTumor)
  .resolve(async (props) => {
    return await python(
      'components.service.kf.main',
      { kargs: [props.inputs.gene_info.ensembl?.gene]}, 
    )
  })
  .build()
