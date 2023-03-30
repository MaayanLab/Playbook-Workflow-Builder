import { MetaNode } from '@/spec/metanode'
import { GeneInfo, GeneInfoFromGeneTerm } from '@/components/service/mygeneinfo'
import { GeneExpressionInTumor } from '@/components/core/input/tumor'
import python from '@/utils/python'
import { GeneTerm } from '@/components/core/input/term'

export const KFTumorExpression = MetaNode('KFTumorExpression')
  .meta({
    label: 'Query KF Gene Expression in Tumor',
    description: 'Use KF API to obtain tumors expressing the given gene',
    pagerank: 1,
  })
  .inputs({ gene_info: GeneInfo })
  .output(GeneExpressionInTumor)
  .resolve(async (props) => {
    return await python(
      'components.service.kf.main',
      { kargs: [props.inputs.gene_info.entrezgene]}, 
    )
  })
  .build()

export const KFTumorExpressionFromGene = MetaNode('KFTumorExpressionFromGene')
  .meta(KFTumorExpression.meta)
  .inputs({ gene: GeneTerm })
  .output(KFTumorExpression.output)
  .resolve(async (props) => {
    const gene_info = await GeneInfoFromGeneTerm.resolve(props)
    return await KFTumorExpression.resolve({ inputs: { gene_info } })
  })
  .build()
