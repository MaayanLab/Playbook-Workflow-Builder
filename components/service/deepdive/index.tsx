import IFrame from "@/app/components/IFrame"
import { GeneTerm } from "@/components/core/term"
import { MetaNode } from "@/spec/metanode"
import { z } from 'zod'

const deepdive_url = 'https://maayanlab.cloud/deepdive'

export const GeneDeepDiveSummary = MetaNode(`GeneDeepDiveSummary`)
  .meta({
    label: 'LLM DeepDive Gene Summary',
    description: 'A LLM DeepDive on a gene',
    external: true,
  })
  .codec(z.string())
  .view(props => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 800 }}>
      <IFrame
        className="flex-grow border-0"
        src={`${deepdive_url}/${props}/`}
      />
    </div>
  ))
  .build()

export const GeneDeepDive = MetaNode(`GeneDeepDive`)
  .meta({
    label: 'LLM DeepDive Gene Summary',
    description: 'Summarize the top 50 most cited abstracts mentioning the gene',
    external: true,
  })
  .inputs({ gene: GeneTerm })
  .output(GeneDeepDiveSummary)
  .resolve(async (props) => {
    return props.inputs.gene
  })
  .story(props => ({
    abstract: `The top 50 most cited paper abstracts that mention ${props.inputs?.gene ?? 'the gene'} were summarized using DeepDive\\ref{DeepDive, https://maayanlab.cloud/deepdive/}.`,
    legend: `A summary of the top 50 most cited paper abstracts that mention ${props.inputs?.gene ?? 'the gene'} using gpt4o-mini\\ref{DeepDive, https://maayanlab.cloud/deepdive/}.`,
  }))
  .build()
