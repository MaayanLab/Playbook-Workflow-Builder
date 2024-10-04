import { MetaNode } from '@/spec/metanode'
import { GeneTerm} from '@/components/core/term'
import { GeneAssociations_HG38 } from '../hg38GeneAssociations'

export const GetGeneFromAssociations = MetaNode('GetGeneFromAssociations')
  .meta({
    label: 'Identify a gene with specified variant ',
    description: 'Use gene associations with variant to identify variants within gene'
  })
  .inputs({ associated: GeneAssociations_HG38})
  .output(GeneTerm)
  .resolve(async (props) => {
        const filtered = props.inputs.associated.filter(geneAssoc => 
            geneAssoc.associations.some(association => 
                association.distance_to_feature === "within gene"
            )
        )
        .map(geneAssoc => geneAssoc.geneId);

        return filtered[0];
    }).story(props => ({
    abstract: `Gene(s) containing the variant(s) from associations are returned.`,
  })).build()