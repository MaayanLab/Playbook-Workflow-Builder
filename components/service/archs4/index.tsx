import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { ScoredGenes, ScoredTissues } from '@/components/core/scored'
import { archs4_icon, weighted_icon } from '@/icons'
import { GeneTerm } from '@/components/core/term'
import read_csv from '@/utils/csv'
import { z } from 'zod'
import { Cell, Column, Table } from '@/app/components/Table'
import { downloadBlob } from '@/utils/download'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import python from '@/utils/python'

async function archs4_tissue_expression({ search, species = 'human', type = 'tissue' }: { search: string, species?: string, type?: string }) {
  const params = new URLSearchParams()
  params.append('search', search)
  params.append('species', species)
  params.append('type', type)
  const req = await fetch(`https://maayanlab.cloud/archs4/search/loadExpressionTissue.php?${params.toString()}`)
  const { values } = read_csv<'id' | 'min' | 'q1' | 'median' | 'q3' | 'max' | 'color'>(await req.text())
  const results = values
    .map(({ id, median }) => ({
      term: id.split('.')[3],
      zscore: (Number(median) - 4.336)/3.921, // these mu/std coefficients were obtained by querying random genes
    }))
    .filter(({ term }) => !!term)
  results.sort((a, b) => b.zscore - a.zscore)
  return results
}

export const ARCHS4TissueExpression = MetaNode('ARCHS4TissueExpression')
  .meta({
    label: 'Query ARCHS4 Median Tissue Expression',
    description: 'Use ARCHS4 API to obtain median tissue expression for the given gene',
    icon: [archs4_icon],
  })
  .inputs({ gene: GeneTerm })
  .output(ScoredTissues)
  .resolve(async (props) => {
    return await archs4_tissue_expression({ search: props.inputs.gene })
  })
  .story(props => ({
    abstract: `Median expression of ${props.inputs?.gene ? props.inputs.gene : 'the gene'} was obtained from ARCHS4\\ref{doi:10.1038/s41467-018-03751-6}.`,
    introduction: `All RNA-seq and ChIP-seq sample and signature search (ARCHS4)\\ref{doi:10.1038/s41467-018-03751-6} is a resource providing access to gene counts uniformly processed from all human and mouse RNA-seq experiments from the Gene Expression Omnibus (GEO) and the Sequence Read Archive (SRA).`,
    methods: `Median expression of the gene in various tissues across experiments from GEO are resolved from the ARCHS4\\ref{doi:10.1038/s41467-018-03751-6} website's REST API. 1000 random genes were queried against the ARCHS4 to get a mean and standard deviation for the median value. The medians returned by the ARCHS4 for each tissue are assigned a z-score based on the pre-computed mean and standard deviations.`,
    legend: `A table showing the median expression of the gene in various tissues from ARCHS4\\ref{doi:10.1038/s41467-018-03751-6}.`,
  }))
  .build()

export const ARCHS4SignatureResults = MetaNode(`ARCHS4SignatureResults`)
  .meta({
    label: 'ARCHS4 Signature Search Results',
    description: 'ARCHS4 Signatures Query Results',
    icon: [archs4_icon, weighted_icon],
  })
  .codec(z.object({ samples: z.array(z.string()) }))
  .view(results =>
    <Table
      height={500}
      cellRendererDependencies={[results]}
      numRows={results.samples.length}
      enableGhostCells
      enableFocusedCell
      downloads={{
        JSON: () => downloadBlob(new Blob([JSON.stringify(results)], { type: 'application/json;charset=utf-8' }), 'data.json'),
        CSV: () => downloadBlob(new Blob([
          [
            `Samples`,
            ...results.samples
          ].join('\n')
        ], { type: 'text/csv;charset=utf-8' }), 'data.csv'),
      }}
    >
      <Column
        name="Samples"
        cellRenderer={row => <Cell key={row+''}><a href={`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${results.samples[row]}`}>{results.samples[row]}</a></Cell>}
      />
    </Table>
  )
  .build()

export const SomeARCHS4SignatureResults = MetaNode(`SomeARCHS4SignatureResults`)
  .meta({
    label: `Select Some GEO Samples`,
    description: `Select some GEO Samples`,
  })
  .codec(z.record(z.string(), z.literal(true)))
  .inputs({ results: ARCHS4SignatureResults })
  .output(ARCHS4SignatureResults)
  .prompt(props => {
    const set = props.inputs.results.samples
    const [selected, setSelected] = React.useState(props.data ? props.data : {} as Record<string, true>)
    React.useEffect(() => {
      if (props.output !== undefined) {
        const selected_ = {} as Record<string, true>
        set.forEach(item => {
          if (!props.output?.samples.includes(item)) {
            selected_[item] = true
          }
        })
        setSelected(selected_)
      }
    }, [props.output, set])
    return (
      <>
        <Table
          height={500}
          cellRendererDependencies={[set]}
          rowHeaderCellRenderer={(row) =>
            <div
              className="text-center block"
              onClick={evt => {
                setSelected(({ [set[row]]: current, ...selected }) => !current ? ({ ...selected, [set[row]]: true }) : selected)
              }}
            >
              <input type="checkbox" checked={!selected[set[row]]} />
            </div>
          }
          numRows={set.length}
          shape={[set.length - Object.keys(selected).length]}
          enableGhostCells
        >
          <Column
            name={"Samples"}
            cellRenderer={row => <Cell key={row+''}>{set[row]}</Cell>}
          />
        </Table>
        <button className="bp5-button bp5-large" onClick={async () => {
          props.submit(selected)
        }}>Submit</button>
      </>
    )
  })
  .resolve(async (props) => {
    const set = new Set()
    props.inputs.results.samples.forEach(item => set.add(item))
    if (Object.keys(props.data).some(key => !set.has(key))) {
      throw new Error('Please select genes from the table')
    }
    return {
      samples: props.inputs.results.samples.filter((item) => !props.data[item]),
    }
  })
  .story(props =>  ({
    abstract: `Some GEO samples were selected for further investigation.`,
    legend: `Some GEO samples taken from ${props.input_refs?.results}.`,
  }))
  .build()

export const ARCHS4SignatureResolve = MetaNode(`ARCHS4SignatureResolve`)
  .meta({
    label: `Fetch GEO Sample Expression`,
    description: `Get the gene count expression data from ARCHS4`,
  })
  .inputs({ results: ARCHS4SignatureResults })
  .output(GeneCountMatrix)
  .resolve(async (props) => await python(
    'components.service.archs4.fetch_samples',
    { kargs: [props.inputs.results.samples] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `Gene expression for published samples was obtained from ARCHS4\\ref{doi:10.1038/s41467-018-03751-6}.`,
    introduction: `All RNA-seq and ChIP-seq sample and signature search (ARCHS4)\\ref{doi:10.1038/s41467-018-03751-6} is a resource providing access to gene counts uniformly processed from all human and mouse RNA-seq experiments from the Gene Expression Omnibus (GEO) and the Sequence Read Archive (SRA).`,
    legend: `A table showing the basic structure and shape of the gene count matrix. Rows represent columns, columns represent genes, and values show the number of mapped reads.`,
  }))
  .build()

export const ARCHS4SignatureSearchT = [
  'Human', 'Mouse',
].map(species => MetaNode(`ARCHS4SignatureSearch[${species}]`)
  .meta({
    label: `ARCHS4 Signature Search`,
    description: `Query ARCHS4 ${species} Signatures`,
    icon: [archs4_icon],
  })
  .inputs({ genes: ScoredGenes })
  .output(ARCHS4SignatureResults)
  .resolve(async (props) => {
    const up: Record<string, true> = {}
    const down: Record<string, true> = {}
    props.inputs.genes.forEach(({ term, zscore }) => {
      if ((typeof zscore === 'number' && zscore >= 1.645) || zscore === 'inf') up[term] = true
      else if ((typeof zscore === 'number' && zscore <= -1.645) || zscore === '-inf') down[term] = true
    })
    const req = await fetch(`https://maayanlab.cloud/rookpy/signature`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Accept': 'application/json',
      },
      body: JSON.stringify({
        type: 'geneset',
        signatureName: 'playbook-workflow-builder.cloud input',
        species: species.toLocaleLowerCase(),
        upgenes: Object.keys(up),
        downgenes: Object.keys(down),
      }),
    })
    if (!req.ok) throw new Error(`ARCHS4 Signature Search Error: ${await req.text()}`)
    return {
      samples: z.object({
        name: z.string(),
        samples: z.array(z.number()),
      }).parse(await req.json()).samples.map(sample => `GSM${sample}`)
    }
  })
  .story(props => ({
    abstract: `Reversers and mimickers from GEO signatures were identified using ARCHS4\\ref{doi:10.1038/s41467-018-03751-6}.`,
    introduction: `All RNA-seq and ChIP-seq sample and signature search (ARCHS4)\\ref{doi:10.1038/s41467-018-03751-6} is a resource providing access to gene counts uniformly processed from all human and mouse RNA-seq experiments from the Gene Expression Omnibus (GEO) and the Sequence Read Archive (SRA).`,
    methods: `Reverser and mimicker signatures from GEO were identified using the ARCHS4\\ref{doi:10.1038/s41467-018-03751-6} website's REST API.`,
    legend: `A table of GEO sample ids for matching signatures identified by ARCHS4\\ref{doi:10.1038/s41467-018-03751-6}.`,
  }))
  .build()
)
