import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { Table, Cell, Column, EditableCell } from '@/app/components/Table'
import { DrugSet } from '@/components/core/input/set'
import { FileC, FilePrompt, FileURL } from '@/components/core/file'
import { downloadBlob } from '@/utils/download'
import { file_icon, file_transfer_icon, gmt_icon } from '@/icons'
import dynamic from 'next/dynamic'
import python from '@/utils/python'
import SafeRender from '@/utils/saferender'

const Bp5Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))

export const DMT = MetaNode(`DMT`)
  .meta({
    label: `Drug Matrix Transpose`,
    description: 'Terms mapped to drugs',
    icon: [gmt_icon],
  })
  .codec(z.record(z.string(), z.object({ description: z.string().optional(), set: z.array(z.string()) })))
  .view(dmt => {
    const items = dict.items(dmt)
    return (
      <Table
        height={500}
        cellRendererDependencies={[items]}
        numRows={items.length}
        enableGhostCells
        enableFocusedCell
        downloads={{
          JSON: () => downloadBlob(new Blob([JSON.stringify(dmt)], { type: 'application/json;charset=utf-8' }), 'data.json'),
          DMT: () => downloadBlob(new Blob([items.map(({ key: term, value: { description, set } }) => [term, description||'', ...set].join('\t')).join('\n')], { type: 'text/tab-separated-values;charset=utf-8' }), 'data.dmt'),
        }}
      >
        <Column
          name="Term"
          cellRenderer={row => <Cell key={row+''}>{items[row].key.toString()}</Cell>}
        />
        <Column
          name="Description"
          cellRenderer={row => <Cell key={row+''}>{items[row].value.description||''}</Cell>}
        />
        <Column
          name="Drug Set"
          cellRenderer={row => <Cell key={row+''}>{items[row].value.set.join('\t')}</Cell>}
        />
      </Table>
    )
  })
  .build()

  export const DMTFromFile = MetaNode('DMTFromFile')
  .meta({
    label: 'Resolve A Drug Matrix Transpose from a File',
    description: 'Ensure a file contains a drug matrix transpose, load it into a standard format',
    icon: [file_transfer_icon],
  })
  .inputs({ file: FileURL })
  .output(DMT)
  .resolve(async (props) => await python(
    'components.data.drug_matrix_transpose.load_drug_matrix_transpose',
    { kargs: [props.inputs.file] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => `The file${props.inputs && props.inputs.file.description ? ` containing ${props.inputs.file.description}` : ''} was loaded as a drug matrix transpose.`)
  .build()

export const DMTFileUpload = MetaNode('DMTFileUpload')
  .meta({
    label: 'Upload A Drug Matrix Transpose',
    description: 'A file containing labeled drug sets',
    tags: {
      Type: {
        File: 1,
        Drug: 1,
      },
      Cardinality: {
        Matrix: 1,
      },
    },
    icon: [file_icon],
  })
  .codec(FileC)
  .inputs()
  .output(DMT)
  .prompt(props => <><FilePrompt {...props} />{props.output ? <SafeRender component={DMT.view} props={props.output} /> : null}</>)
  .resolve(async (props) => await python(
    'components.data.drug_matrix_transpose.load_drug_matrix_transpose',
    { kargs: [props.data] },
    message => props.notify({ type: 'info', message }),
  ))
  .story(props => `A drug matrix transpose${props.data && props.data.description ? ` containing ${props.data.description}` : ''} was uploaded.`)
  .build()

export const DMTUnion = MetaNode('DMTUnion')
  .meta({
    label: `Compute Union Drug Set`,
    description: 'Find the union set of all drugs in the DMT'
  })
  .inputs({ dmt: DMT })
  .output(DrugSet)
  .resolve(async (props) => {
    return { set: array.unique(dict.values(props.inputs.dmt).flatMap(({ set: drugset }) => drugset)) }
  })
  .story(props =>
    `All the identified drug sets were combined using the union set operation.`
  )
  .build()

export const DMTIntersection = MetaNode('DMTIntersection')
  .meta({
    label: `Compute Intersection Drug set`,
    description: 'Find the intersecting set of all drugs in the DMT'
  })
  .inputs({ dmt: DMT })
  .output(DrugSet)
  .resolve(async (props) => {
    return dict.values(props.inputs.dmt).reduce(({ set: A }, { set: B }) => ({ set: array.intersection(A, B) }))
  })
  .story(props =>
    `A consensus drug set was created using the set intersection operation.`
  )
  .build()

export const DMTConsensus = MetaNode('DMTConsensus')
  .meta({
    label: `Compute Consensus Drug Set`,
    description: 'Find drugs which appear in more than one set'
  })
  .inputs({ dmt: DMT })
  .output(DrugSet)
  .resolve(async (props) => {
    const counts: Record<string, number> = {}
    dict.values(props.inputs.dmt)
      .forEach(({ set }) =>
        set.forEach(drug =>
          counts[drug] = (counts[drug]||0)+1
        )
      )
    return {
      set: dict.items(counts)
        .filter(({ value }) => value > 1)
        .map(({ key }) => key as string)
    }
  })
  .story(props =>
    `A consensus drug set was created by only retaining drugs that appear in at least two sets.`
  )
  .build()

export const DrugSetsToDMT = MetaNode('DrugSetsToDMT')
  .meta({
    label: `Assemble DMT from Drug Sets`,
    description: 'Group multiple independently generated drug sets into a single DMT'
  })
  .codec(z.object({
    terms: z.record(z.union([z.string(), z.number()]).transform(key => +key), z.string()),
    descriptions: z.record(z.union([z.string(), z.number()]).transform(key => +key), z.string()),
  }))
  .inputs({ sets: [DrugSet] })
  .output(DMT)
  .prompt(props => {
    const [terms, setTerms] = React.useState(() => dict.init(array.arange(props.inputs.sets.length).map(key => ({ key, value: props.inputs.sets[key].description||'' }))))
    const [descriptions, setDescriptions] = React.useState(() => dict.init(array.arange(props.inputs.sets.length).map(key => ({ key, value: '' }))))
    React.useEffect(() => {
      if (!props.data) return
      setTerms(props.data.terms ? props.data.terms : {} as Record<number, string>)
      setDescriptions(props.data.descriptions ? props.data.descriptions : {} as Record<number, string>)
    }, [props.data])
    return (
      <div>
        <Table
          height={500}
          cellRendererDependencies={[props.inputs.sets, terms, descriptions]}
          numRows={props.inputs.sets.length}
          enableGhostCells
        >
          <Column
            name="Term"
            cellRenderer={row => <EditableCell
              key={row+''}
              value={terms[row]}
              editableTextProps={{placeholder: `Drug set from path ${row}`}}
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
            name="Drug Set"
            cellRenderer={row => <Cell key={row+''}>{props.inputs.sets[row].set.join('\t')}</Cell>}
          />
        </Table>
        <Bp5Button
          large
          type="submit"
          text="Submit"
          rightIcon="bring-data"
          onClick={() => props.submit({ terms, descriptions })}
        />
      </div>
    )
  })
  .resolve(async (props) => {
    const { terms, descriptions } = props.data
    if (props.inputs.sets.length !== Object.keys(terms).length) {
      throw new Error('Please confirm the terms')
    }
    return dict.init(
      array.arange(props.inputs.sets.length)
        .map(i => ({
          key: terms[i],
          value: {
            description: descriptions[i],
            set: props.inputs.sets[i].set,
          }
        }))
    )
  })
  .story(props =>
    `The drug sets collected were combined into one drug set library.`
  )
  .build()

export const DMTConcatenate = MetaNode('DMTConcatenate')
  .meta({
    label: `Concatenate DMTs`,
    description: 'Join several DMTs into one'
  })
  .inputs({ dmts: [DMT] })
  .output(DMT)
  .resolve(async (props) => {
    return dict.init(props.inputs.dmts.flatMap(dict.items))
  })
  .story(props =>
    `The identified drug sets were combined into one drug set library.`
  )
  .build()
