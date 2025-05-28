import { ScoredDrugs } from "@/components/core/scored"
import { DrugSet } from "@/components/core/set"
import { MetaNode } from "@/spec/metanode"
import python from "@/utils/python"

export const FilterFDAApprovedSetDrugs = MetaNode(`FilterFDAApprovedDrugs[Set[Drug]]`)
  .meta({
    label: 'Extract FDA Approved Drugs',
    description: 'Filter the drug set with PubChem APIs',
    external: true,
  })
  .inputs({ drugs: DrugSet })
  .output(DrugSet)
  .resolve(async (props) => await python(
      'components.service.pubchem.filter_fda_set',
      { kargs: [props.inputs.drugs] },
      message => props.notify({ type: 'info', message }),
    )
  )
  .story(props => ({
    abstract: `The drug set was filtered by FDA Approved Drugs with the help of PubChem APIs\\ref{doi:10.1093/nar/gkac956}.`,
    introduction: `PubChem is an online database of chemical information which can be queried via REST API\\ref{doi:10.1093/nar/gkac956}.`,
    methods: `Drugs in ${props.input_refs?.drugs} were queried against PubChem's API\\ref{doi:10.1093/nar/gkac956} and filtered by those which could be located and determined to be FDA approved to produce ${props.output_ref}.`,
    tableLegend: `A table showing the set of drugs from ${props.input_refs?.drugs} which are FDA approved.`,
  }))
  .build()

  export const FilterFDAApprovedScoredDrugs = MetaNode(`FilterFDAApprovedDrugs[Scored[Drug]]`)
  .meta({
    label: 'Extract FDA Approved Drugs',
    description: 'Filter the drugs with PubChem APIs',
    external: true,
  })
  .inputs({ drugs: ScoredDrugs })
  .output(ScoredDrugs)
  .resolve(async (props) => await python(
      'components.service.pubchem.filter_fda_scored',
      { kargs: [props.inputs.drugs] },
      message => props.notify({ type: 'info', message }),
    )
  )
  .story(props => ({
    abstract: `The drugs were filtered by FDA Approved Drugs with the help of PubChem APIs\\ref{doi:10.1093/nar/gkac956}.`,
    introduction: `PubChem is an online database of chemical information which can be queried via REST API\\ref{doi:10.1093/nar/gkac956}.`,
    methods: `Drugs in ${props.input_refs?.drugs} were queried against PubChem's API\\ref{doi:10.1093/nar/gkac956} and filtered by those which could be located and determined to be FDA approved to produce ${props.output_ref}.`,
    tableLegend: `A table showing the set of drugs from ${props.input_refs?.drugs} which are FDA approved.`,
  }))
  .build()
