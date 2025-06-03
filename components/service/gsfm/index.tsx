import IFrame from "@/app/components/IFrame"
import { gsfm_icon } from '@/icons'
import { GeneTerm } from "@/components/core/term"
import { MetaNode } from "@/spec/metanode"
import { z } from 'zod'

const gsfm_url = 'https://gsfm.maayanlab.cloud'

export const GSFMSummary = MetaNode(`GSFMSummary`)
  .meta({
    label: 'GSFM Predictions',
    description: 'Predictions by Gene Set Foundation Model (GSFM)\\ref{doi:10.1101/2025.05.30.657124} trained on millions of gene sets from literature',
    icon: [gsfm_icon],
    external: true,
  })
  .codec(z.string())
  .view(props => (
    <div className="flex-grow flex flex-row m-0" style={{ minHeight: 800 }}>
      <IFrame
        className="flex-grow border-0"
        src={`${gsfm_url}/gene/${props}?embed`}
      />
    </div>
  ))
  .build()

export const GSFM = MetaNode(`GSFM`)
  .meta({
    label: 'GSFM Predictions',
    description: 'Obtain gene function predictions from Gene Set Foundation Model (GSFM)',
    icon: [gsfm_icon],
    external: true,
  })
  .inputs({ gene: GeneTerm })
  .output(GSFMSummary)
  .resolve(async (props) => {
    return props.inputs.gene
  })
  .story(props => ({
    abstract: `Gene function predictions for ${props.inputs?.gene ?? 'the gene'} were made using GSFM\\ref{doi:10.1101/2025.05.30.657124}.`,
    legend: `Gene function predictions for ${props.inputs?.gene ?? 'the gene'} by GSFM\\ref{doi:10.1101/2025.05.30.657124}.`,
  }))
  .build()
