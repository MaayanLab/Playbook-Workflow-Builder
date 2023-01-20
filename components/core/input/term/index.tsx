import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { Gene, Drug, Primative, Metabolite, Pathway, Phenotype, Tissue, Disease } from '@/components/core/input/primitives'
import dynamic from 'next/dynamic'
import { input_icon } from '@/icons'

const Suggest2 = dynamic(() => import('@blueprintjs/select').then(({ Suggest2 }) => Suggest2))
const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))
const MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))

const Term_T = (T: Primative) => MetaNode.createData(`Term[${T.name}]`)
  .meta({
    label: T.label,
    description: `${T.label} Term`,
    color: T.color,
    icon: T.icon,
    ...(T.extra?.term?.meta || {}),
  })
  .codec(z.string())
  .view(term => {
    return <div>{T.label}: {term}</div>
  })
  .build()

export const DiseaseTerm = Term_T(Disease)
export const DrugTerm = Term_T(Drug)
export const GeneTerm = Term_T(Gene)
export const MetaboliteTerm = Term_T(Metabolite)
export const PathwayTerm = Term_T(Pathway)
export const PhenotypeTerm = Term_T(Phenotype)
export const TissueTerm = Term_T(Tissue)

const itemRenderer = (item: unknown, { modifiers: { active, disabled }, handleClick }: { modifiers: { active: boolean, disabled: boolean }, handleClick: React.MouseEventHandler }) => (
  <MenuItem
    key={item+''}
    text={item+''}
    onClick={handleClick}
    active={active}
    disabled={disabled}
  />
)
const createNewItemRenderer = (item: string, active: boolean, handleClick: React.MouseEventHandler<HTMLElement>) => (
  <MenuItem
    key={item}
    text={item}
    onClick={handleClick}
    active={active}
  />
)
const createNewItemFromQuery = (item: unknown) => item+''
const inputValueRenderer = (item: unknown) => item+''

const Input_Term_T = (T: Primative, Term_T: typeof GeneTerm) => MetaNode.createProcess(`Input[${T.name}]`)
  .meta({
    label: `${T.label} Input`,
    description: `Start with a ${T.label} term`,
    icon: [input_icon],
    tags: {
      Type: {
        [T.label]: 1,
      },
      Cardinality: {
        Term: 1,
      },
    },
    ...(T.extra?.term?.meta || {}),
  })
  .inputs()
  .output(Term_T)
  .prompt(props => {
    const [item, setItem] = React.useState('')
    const [query, setQuery] = React.useState('')
    const { items, error } = T.extra?.term?.autocomplete !== undefined ? T.extra.term.autocomplete(query) : { items: [], error: undefined }
    if (error) console.warn(error)
    React.useEffect(() => { setItem(props.output || '') }, [props.output])
    return (
      <div>
        <Suggest2
          fill
          closeOnSelect
          selectedItem={item}
          createNewItemFromQuery={createNewItemFromQuery}
          onItemSelect={item => setItem(item as string)}
          inputValueRenderer={inputValueRenderer}
          itemRenderer={itemRenderer}
          createNewItemRenderer={createNewItemRenderer}
          noResults={
            <MenuItem
              key="No results"
              text="No results"
              disabled
            />
          }
          items={items}
          inputProps={{ leftIcon: 'search', placeholder: `Search ${T.label}...` }}
          popoverProps={{ minimal: true }}
          onQueryChange={q => setQuery(q)}
        />
        <Button
          large
          type="submit"
          text="Submit"
          rightIcon="send-to-graph"
          onClick={evt => props.submit(item)}
        />
        {T.extra?.term?.meta?.example !== undefined ?
          <Button
            large
            text="Example"
            rightIcon="bring-data"
            onClick={evt => {
              if (T.extra?.term?.meta?.example !== undefined) {
                props.submit(T.extra.term.meta.example)
              }
            }}
          />
          : null}
      </div>
    )
  })
  .build()

export const InputDiseaseTerm = Input_Term_T(Disease, DiseaseTerm)
export const InputDrugTerm = Input_Term_T(Drug, DrugTerm)
export const InputGeneTerm = Input_Term_T(Gene, GeneTerm)
export const InputMetaboliteTerm = Input_Term_T(Metabolite, MetaboliteTerm)
export const InputPathwayTerm = Input_Term_T(Pathway, PathwayTerm)
export const InputPhenotypeTerm = Input_Term_T(Phenotype, PhenotypeTerm)
export const InputTissueTerm = Input_Term_T(Tissue, TissueTerm)
