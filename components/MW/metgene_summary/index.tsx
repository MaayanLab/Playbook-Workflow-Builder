import { z } from 'zod'
import { MetaNode } from '@/spec/metanode'
import { metgene_icon, plot_icon } from '@/icons'
import { GeneTerm } from '@/components/core/input/term'

export const MetGeneSummary = MetaNode('MetGeneSummary')
  .meta({
    label: `MetGENE Summary`,
    description: 'A dashboard for reviewing gene-centric information for a given gene from metabolomics',
    icon: [metgene_icon, plot_icon],
  })
  .codec(z.object({ gene: z.string() }))
  .view(value => {
    return (
      <div className="flex-grow flex flex-row m-0" style={{ minHeight: 500 }}>
        <iframe
          className="flex-grow border-0"
          src={`https://bdcw.org/MetGENE/summary.php?species=hsa&GeneIDType=SYMBOL&anatomy=NA&disease=NA&phenotype=NA&GeneInfoStr=${encodeURIComponent(value.gene)}`}
        />
      </div>
    )
  })
  .build()

export const MetGeneSearch = MetaNode('MetGeneSearch')
  .meta({
    label: `MetGENE Search`,
    description: 'Identify gene-centric information from Metabolomics.',
    icon: [metgene_icon],
    pagerank: 2,
  })
  .inputs({ gene: GeneTerm })
  .output(MetGeneSummary)
  .resolve(async (props) => {
    return { gene: props.inputs.gene }
  })
  .build()
