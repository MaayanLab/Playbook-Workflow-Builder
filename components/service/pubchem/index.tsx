import { ScoredDrugs } from "@/components/core/scored"
import { DrugSet } from "@/components/core/set"
import { MetaNode } from "@/spec/metanode"
import python from "@/utils/python"

export const FilterFDAApprovedSetDrugs = MetaNode(`FilterFDAApprovedDrugs[Set[Drug]]`)
  .meta({
    label: 'Extract FDA Approved Drugs',
    description: 'Filter the drug set with PubChem APIs',
  })
  .inputs({ drugs: DrugSet })
  .output(DrugSet)
  .resolve(async (props) => await python(
      'components.service.pubchem.filter_fda_set',
      { kargs: [props.inputs.drugs] },
      message => props.notify({ type: 'info', message }),
    )
  )
  .story(props => `The drug set was filtered by FDA Approved Drugs with the help of PubChem APIs [\\ref{doi:10.1093/nar/gkac956}].`)
  .build()

  export const FilterFDAApprovedScoredDrugs = MetaNode(`FilterFDAApprovedDrugs[Scored[Drug]]`)
  .meta({
    label: 'Extract FDA Approved Drugs',
    description: 'Filter the drugs with PubChem APIs',
  })
  .inputs({ drugs: ScoredDrugs })
  .output(ScoredDrugs)
  .resolve(async (props) => await python(
      'components.service.pubchem.filter_fda_scored',
      { kargs: [props.inputs.drugs] },
      message => props.notify({ type: 'info', message }),
    )
  )
  .story(props => `The drugs were filtered by FDA Approved Drugs with the help of PubChem APIs [\\ref{doi:10.1093/nar/gkac956}].`)
  .build()
