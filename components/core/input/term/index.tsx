import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { z } from 'zod'
import { Gene, Drug, Primative } from '@/components/core/input/primitives'
import dynamic from 'next/dynamic'

const Suggest2 = dynamic(() => import('@blueprintjs/select').then(({ Suggest2 }) => Suggest2))
const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))
const MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))

const Term_T = (T: Primative) => MetaNode.createData(`Term[${T.name}]`)
  .meta({
    label: T.label,
    description: `${T.label} Term`,
    color: T.color,
    example: T.examples.term,
  })
  .codec(z.string())
  .view(term => {
    return <div>{T.label}: {term}</div>
  })
  .build()

export const GeneTerm = Term_T(Gene)
export const DrugTerm = Term_T(Drug)


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
  })
  .inputs()
  .output(Term_T)
  .prompt(props => {
    const [item, setItem] = React.useState('')
    const [query, setQuery] = React.useState('')
    const items: string[] = []
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
        <div>
          <Button
            large
            type="submit"
            text="Submit"
            rightIcon="send-to-graph"
            onClick={evt => props.submit(item)}
          />
        </div>
        <div>
          <Button
            large
            text="Example"
            rightIcon="bring-data"
            onClick={evt => props.submit(T.examples.term)}
          />
        </div>
      </div>
    )
  })
  .build()

export const InputGeneTerm = Input_Term_T(Gene, GeneTerm)
export const InputDrugTerm = Input_Term_T(Drug, DrugTerm)
