import { z } from 'zod'
import { MetaNode } from '@/spec/metanode'
import { metgene_icon, plot_icon } from '@/icons'
import { GeneTerm } from '@/components/core/term'
import dynamic from 'next/dynamic'

const IFrame = dynamic(() => import('@/app/components/IFrame'))

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
        <IFrame
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
  .story(props => ({
    abstract: `The gene was searched with the MetGENE tool providing pathways, reactions, metabolites, and studies from the Metabolomics Workbench\\ref{MetGENE, https://sc-cfdewebdev.sdsc.edu/MetGENE/metGene.php}.`,
    introduction: `MetGENE is a information retrieval tool that connects a gene or a set of genes to metabolomic studies in the Metabolomic Workbench. It uses a knowledge based approach where the gene is connected to pathways it regulates, followed by reactions within the pathways and metabolites particicpating in the reactions. The metabolites are connected to studies in Metabolomics Workbench.`,
    methods: `Given a gene or a geneset, MetGENE provides a summary of pathways inw hich the gene participates, reactions corresponding to those pathways, metabolites participating in those reactions and the metabolomic studies in which each of those metabolites are measured. It provides a REST API for obtaining the summary is\\ref{https://bdcw.org/MetGENE/rest/summary/species/hsa/GeneIDType/SYMBOL/GeneInfoStr/HK1/anatomy/NA/disease/NA/phenotype/NA/viewType/json} and the output is a table comprising of the number of pathways, reactions, metabolites and studies corresponding to the gene.`,
    legend: `The list of results summarizes the number of pathways in which a gene or gene set participates, as well as the total number of unique reactions in those pathways, metabolites involved in the reactions, and studies in which those metabolites are present.`,
  }))
  .build()
