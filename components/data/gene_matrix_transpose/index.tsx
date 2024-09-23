import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { Table, Cell, Column, EditableCell } from '@/app/components/Table'
import { GeneSet } from '@/components/core/set'
import { FileC, FilePrompt, FileURL } from '@/components/core/file'
import { downloadBlob } from '@/utils/download'
import { file_icon, file_transfer_icon, gmt_icon } from '@/icons'
import dynamic from 'next/dynamic'
import python from '@/utils/python'
import SafeRender from '@/utils/saferender'
import { clientLoadExample } from '@/components/data/gene_matrix_transpose/api/Aging_Perturbations_from_GEO_Brain_up.gmt/client'

const Bp5Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export const GMT = MetaNode(`GMT`)
  .meta({
    label: `Gene Matrix Transpose`,
    description: 'Terms mapped to genes',
    icon: [gmt_icon],
  })
  .codec(z.record(z.string(), z.object({ description: z.string().optional(), set: z.array(z.string()) })))
  .view(gmt => {
    const gmt_items = dict.items(gmt)
    return (
      <Table
        height={500}
        cellRendererDependencies={[gmt_items]}
        numRows={gmt_items.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(gmt)], { type: 'application/json;charset=utf-8' }), 'data.json'),
          GMT: () => downloadBlob(new Blob([gmt_items.map(({ key: term, value: { description, set } }) => [term, description||'', ...set].join('\t')).join('\n')], { type: 'text/tab-separated-values;charset=utf-8' }), 'data.gmt'),
        }}
      >
        <Column
          name="Term"
          cellRenderer={row => <Cell key={row+''}>{gmt_items[row].key.toString()}</Cell>}
        />
        <Column
          name="Description"
          cellRenderer={row => <Cell key={row+''}>{gmt_items[row].value.description||''}</Cell>}
        />
        <Column
          name="Gene Set"
          cellRenderer={row => <Cell key={row+''}>{gmt_items[row].value.set.join('\t')}</Cell>}
        />
      </Table>
    )
  })
  .build()

export const GMTFromFile = MetaNode('GMTFromFile')
  .meta({
    label: 'Resolve A Gene Matrix Transpose from a File',
    description: 'Ensure a file contains a gene matrix transpose, load it into a standard format',
    icon: [file_transfer_icon],
  })
  .inputs({ file: FileURL })
  .output(GMT)
  .resolve(async (props) => await python(
    'components.data.gene_matrix_transpose.load_gene_matrix_transpose',
    { kargs: [props.inputs.file] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''} was loaded as a gene matrix transpose.`,
    introduction: `The gene matrix transpose (GMT) format\\ref{Gene Matrix Transpose file format, https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29} is an efficient sparse matrix format well suited for gene set libraries.`,
    legend: `A collection of gene sets in GMT format.`,
  }))
  .build()

export const GMTFileUpload = MetaNode('GMTFileUpload')
  .meta({
    label: 'Upload A Gene Matrix Transpose',
    description: 'A file containing labeled gene sets',
    tags: {
      Type: {
        File: 1,
        Gene: 1,
      },
      Cardinality: {
        Matrix: 1,
      },
    },
    icon: [file_icon],
  })
  .codec(FileC)
  .inputs()
  .output(GMT)
  .prompt(props => <><FilePrompt {...props} example={clientLoadExample} />{props.output ? <SafeRender component={GMT.view} props={props.output} /> : null}</>)
  .resolve(async (props) => await python(
    'components.data.gene_matrix_transpose.load_gene_matrix_transpose',
    { kargs: [props.data] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => ({
    abstract: `The gene matrix transpose${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`,
    introduction: `The gene matrix transpose (GMT) format\\ref{Gene Matrix Transpose file format, https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29} is an efficient sparse matrix format well suited for gene set libraries.`,
    legend: `A collection of gene sets in GMT format.`,
  }))
  .build()

export const GMTUnion = MetaNode('GMTUnion')
  .meta({
    label: `Compute Union Gene Set`,
    description: 'Find the union set of all genes in the GMT'
  })
  .inputs({ gmt: GMT })
  .output(GeneSet)
  .resolve(async (props) => {
    return { set: array.unique(dict.values(props.inputs.gmt).flatMap(({ set: geneset }) => geneset)) }
  })
  .story(props => ({
    abstract: `All the identified gene sets were combined using the union set operation.`,
    introduction: `The gene matrix transpose (GMT) format\\ref{Gene Matrix Transpose file format, https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29} is an efficient sparse matrix format well suited for gene set libraries.`,
    methods: `All genes present in all gene sets from ${props.input_refs?.gmt} were aggregated into a single gene set.`,
    legend: `The union of all genes across all gene sets in ${props.input_refs?.gmt}.`,
  }))
  .build()

export const GMTIntersection = MetaNode('GMTIntersection')
  .meta({
    label: `Compute Intersection Gene Set`,
    description: 'Find the intersecting set of all genes in the GMT'
  })
  .inputs({ gmt: GMT })
  .output(GeneSet)
  .resolve(async (props) => {
    return dict.values(props.inputs.gmt).reduce(({ set: A }, { set: B }) => ({ set: array.intersection(A, B) }))
  })
  .story(props => ({
    abstract: `A consensus gene set was created using the set intersection operation.`,
    introduction: `The gene matrix transpose (GMT) format\\ref{Gene Matrix Transpose file format, https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29} is an efficient sparse matrix format well suited for gene set libraries.`,
    methods: `Genes that appear in all gene sets in ${props.input_refs?.gmt} were aggregated into a single gene set.`,
    legend: `The genes that appear in all gene sets in ${props.input_refs?.gmt}.`,
  }))
  .build()

export const GMTConsensus = MetaNode('GMTConsensus')
  .meta({
    label: `Compute Consensus Gene Set`,
    description: 'Find genes which appear in more than one set'
  })
  .inputs({ gmt: GMT })
  .output(GeneSet)
  .resolve(async (props) => {
    const gene_counts: Record<string, number> = {}
    dict.values(props.inputs.gmt)
      .forEach(({ set: geneset }) =>
        geneset.forEach(gene =>
          gene_counts[gene] = (gene_counts[gene]||0)+1
        )
      )
    return {
      set: dict.items(gene_counts)
        .filter(({ value }) => value > 1)
        .map(({ key }) => key as string)
    }
  })
  .story(props => ({
    abstract: `A consensus gene set was created by only retaining genes that appear in at least two sets.`,
    introduction: `The gene matrix transpose (GMT) format\\ref{Gene Matrix Transpose file format, https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29} is an efficient sparse matrix format well suited for gene set libraries.`,
    methods: `Genes that appear in at least 2 of the gene sets in ${props.input_refs?.gmt} were collected into a set.`,
    legend: `The genes that appear in at least 2 of the gene sets in ${props.input_refs?.gmt}.`,
  }))
  .build()

export const GenesetsToGMT = MetaNode('GenesetsToGMT')
  .meta({
    label: `Assemble GMT from Gene Sets`,
    description: 'Group multiple independently generated gene sets into a single GMT'
  })
  .codec(z.object({
    terms: z.record(z.union([z.string(), z.number()]).transform(key => +key), z.string()),
    descriptions: z.record(z.union([z.string(), z.number()]).transform(key => +key), z.string()),
  }))
  .inputs({ genesets: [GeneSet] })
  .output(GMT)
  .prompt(props => {
    const [terms, setTerms] = React.useState(() => dict.init(array.arange(props.inputs.genesets.length).map(key => ({ key, value: props.inputs.genesets[key].description||'' }))))
    const [descriptions, setDescriptions] = React.useState(() => dict.init(array.arange(props.inputs.genesets.length).map(key => ({ key, value: '' }))))
    React.useEffect(() => {
      if (!props.data) return
      setTerms(props.data.terms ? props.data.terms : {} as Record<number, string>)
      setDescriptions(props.data.descriptions ? props.data.descriptions : {} as Record<number, string>)
    }, [props.data])
    return (
      <>
        <Table
          height={500}
          cellRendererDependencies={[props.inputs.genesets, terms, descriptions]}
          numRows={props.inputs.genesets.length}
          enableGhostCells
        >
          <Column
            name="Term"
            cellRenderer={row => <EditableCell
              key={row+''}
              value={terms[row]}
              editableTextProps={{placeholder: `Gene set from path ${row}`}}
              onChange={value => setTerms(terms => ({ ...terms, [row]: value }))}
            />}
          />
          <Column
            name="Description"
            cellRenderer={row => <EditableCell
              key={row+''}
              value={descriptions[row]}
              editableTextProps={{placeholder: `Some optional description`}}
              onChange={value => setDescriptions(descriptions => ({ ...descriptions, [row]: value }))}
            />}
          />
          <Column
            name="Gene Set"
            cellRenderer={row => <Cell key={row+''}>{props.inputs.genesets[row].set.join('\t')}</Cell>}
          />
        </Table>
        <Bp5Button
          large
          type="submit"
          text="Submit"
          rightIcon="cloud-upload"
          onClick={() => props.submit({ terms, descriptions })}
        />
      </>
    )
  })
  .resolve(async (props) => {
    const { terms, descriptions } = props.data
    if (props.inputs.genesets.length !== Object.keys(terms).length) {
      throw new Error('Please confirm the terms')
    }
    return dict.init(
      array.arange(props.inputs.genesets.length)
        .map(i => ({
          key: terms[i],
          value: {
            description: descriptions[i],
            set: props.inputs.genesets[i].set,
          }
        }))
    )
  })
  .story(props => ({
    abstract: `The gene sets collected were combined into one gene set library.`,
    introduction: `The gene matrix transpose (GMT) format\\ref{Gene Matrix Transpose file format, https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29} is an efficient sparse matrix format well suited for gene set libraries.`,
    methods: `Labels were created for each of the gene sets ${props.input_refs?.genesets} to construct a GMT.`,
    legend: `The gene sets from ${props.input_refs?.genesets} collected into a GMT.`,
  }))
  .build()

export const GMTConcatenate = MetaNode('GMTConcatenate')
  .meta({
    label: `Concatenate GMTs`,
    description: 'Join several GMTs into one'
  })
  .inputs({ gmts: [GMT] })
  .output(GMT)
  .resolve(async (props) => {
    return dict.init(props.inputs.gmts.flatMap(dict.items))
  })
  .story(props => ({
    abstract: `Multiple GMTs were combined into one GMT.`,
    introduction: `The gene matrix transpose (GMT) format\\ref{Gene Matrix Transpose file format, https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29} is an efficient sparse matrix format well suited for gene set libraries.`,
    methods: `A joint GMT was constructed by stacking GMTs${props.input_refs?.gmts ? ` from ${(props.input_refs?.gmts as string[]).join(', ')}` : ''}.`,
    legend: `The GMTs${props.input_refs?.gmts ? ` from ${(props.input_refs?.gmts as string[]).join(', ')}` : ''} collected into a joint GMT.`,
     
  }))
  .build()
