import { NextRouter, useRouter } from 'next/router'
import type { Metapath } from '@/app/fragments/graph/types'
import krg from '@/app/krg'
import React from 'react'
import { MetaNode, MetaNodeDataType, MetaNodePromptType, MetaNodeResolveType } from '@/spec/metanode'
import { z } from 'zod'
import { Intent } from '@blueprintjs/core'
import dynamic from 'next/dynamic'
import * as dict from '@/utils/dict'

const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))
const FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const ControlGroup = dynamic(() => import('@blueprintjs/core').then(({ ControlGroup }) => ControlGroup))
const InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))
const TextArea = dynamic(() => import('@blueprintjs/core').then(({ TextArea }) => TextArea))
const MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))
const Suggest2 = dynamic(() => import('@blueprintjs/select').then(({ Suggest2 }) => Suggest2))

const Suggestion = MetaNode.createData('Suggestion')
  .meta({
    label: 'Suggestion',
    description: 'An actual suggestion',
  })
  .codec(z.object({
    name: z.string(),
    inputs: z.string(),
    output: z.string(),
    author_name: z.string(),
    author_email: z.string(),
    author_org: z.string(),
    description: z.string()
  }))
  .view(suggestion => (
    <div>
      <div className="bp4-card">
        <h4 className="bp4-heading">{suggestion.name}</h4>
        <p className="bp4-text-large">{suggestion.description}</p>
        <div className="bp4-callout">
          <h5 className="bp4-heading">Author</h5>
          {suggestion.author_name} &lt;{suggestion.author_email}&gt; ({suggestion.author_org})
        </div>
        {/* <br />
        <h5 className="bp4-heading"><i>Comments</i></h5>
        <div className="bp4-callout">
          <p className="bp4-text-large">My comment is important</p>
          <h5 className="bp4-heading">Author</h5>
          {suggestion.author_name} &lt;{suggestion.author_email}&gt; ({suggestion.author_org})
        </div> */}
      </div>
    </div>
  ))
  .build()

export function SuggestionEdges(input?: MetaNodeDataType) {
  const suggestion_edges: Array<MetaNodePromptType<any> | MetaNodeResolveType<any>> = []
  if (input === undefined) {
    suggestion_edges.push(
      MetaNode.createProcess(`SuggestDataType`)
        .meta({
          label: 'Suggest a core data type',
          description: `This would be usable as an initial or intermediary data`,
        })
        .inputs({})
        .output(Suggestion)
        .prompt((props) => <></>)
        .build()
    )
  }
  if (input !== undefined) {
    suggestion_edges.push(
      MetaNode.createProcess(`SuggestInteractiveEdge[${input.spec}]`)
        .meta({
          label: 'Suggest a visualization method',
          description: `This would visualize the ${input.meta.label || input.spec}. Provide a description about what should be here.`,
        })
        .inputs({ input })
        .output(Suggestion)
        .prompt((props) => <></>)
        .build()
    )
    suggestion_edges.push(
      MetaNode.createProcess(`SuggestResolveEdge[${input.spec}]`)
        .meta({
          label: 'Suggest an algorithm or data transformation method',
          description: `This would transform the ${input.meta.label || input.spec}. Provide a description about what should be here.`,
        })
        .inputs({ input })
        .output(Suggestion)
        .prompt((props) => <></>)
        .build()
    )
  }
  return suggestion_edges.map(element => ({
    ...element,
    onClick: ({ router, id, head }: { router: NextRouter, id: string, head: Metapath }) => {
      router.push(`/graph/${id}/${head ? head.id !== id ? `/node/${head.id}` : '' : '/node/start'}/suggest`)
    }
  }))
}

export default function Suggest({ id, head }: { id: string, head: Metapath }) {
  const router = useRouter()
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  const input = processNode ? processNode.output : undefined
  const [suggestion, setSuggestion] = React.useState({
    name: '',
    inputs: input ? input.spec as string : '',
    output: '',
    author_name: '',
    author_email: '',
    author_org: '',
    description: '',
  })
  return (
    <div>
      <ControlGroup fill vertical>
        <FormGroup
          label="Component Name"
          labelInfo="(required)"
          helperText="A succict name which describes this component"
        >
          <InputGroup
            placeholder="Some component"
            value={suggestion.name}
            onChange={evt => {
              setSuggestion(({ ...suggestion }) => ({ ...suggestion, name: evt.target.value }))
            }}
            leftIcon="label"
          />
        </FormGroup>
        {suggestion.inputs ?
          <ControlGroup fill>
            <FormGroup
              label="Component Inputs"
              labelInfo="(required)"
              helperText="A the inputs to this component"
            >
              <InputGroup
                placeholder="Some component"
                value={suggestion.inputs}
                readOnly
                onChange={evt => {
                  setSuggestion(({ ...suggestion }) => ({ ...suggestion, inputs: evt.target.value }))
                }}
                leftIcon="many-to-one"
              />
            </FormGroup>
            <FormGroup
              label="Component Output"
              labelInfo="(required)"
              helperText="The output of this component"
            >
              <Suggest2
                fill
                closeOnSelect
                selectedItem={suggestion.output}
                inputValueRenderer={item => item+''}
                itemRenderer={(item, { modifiers: { active, disabled }, handleClick }: { modifiers: { active: boolean, disabled: boolean }, handleClick: React.MouseEventHandler }) =>
                  <MenuItem
                    key={item+''}
                    text={item+''}
                    onClick={handleClick}
                    active={active}
                    disabled={disabled}
                  />
                }
                createNewItemFromQuery={(item: string) => item}
                createNewItemRenderer={(item: string, active: boolean, handleClick: React.MouseEventHandler<HTMLElement>) =>
                  <MenuItem
                    key={item}
                    text={item}
                    active={active}
                    onClick={handleClick}
                  />
                }
                inputProps={{ leftIcon: 'flow-end', placeholder: `Search components` }}
                items={krg.getDataNodes().map(node => node.spec)}
                onItemSelect={output => {
                  setSuggestion(({ ...suggestion }) => ({ ...suggestion, output: output+'' }))
                }}
                popoverProps={{ minimal: true }}
              />
            </FormGroup>
          </ControlGroup>
        : null}
      </ControlGroup>
      <FormGroup
        label="Authorship Information"
        labelInfo="(required)"
        helperText="Let us know who you are and how to contact you"
      >
        <ControlGroup fill vertical>
          <ControlGroup fill>
            <InputGroup
              type="text"
              placeholder="Name"
              value={suggestion.author_name}
              onChange={evt => {
                setSuggestion(({ ...suggestion }) => ({ ...suggestion, author_name: evt.target.value }))
              }}
              leftIcon="person"
            />
            <InputGroup
              type="email"
              placeholder="Email"
              value={suggestion.author_email}
              onChange={evt => {
                setSuggestion(({ ...suggestion }) => ({ ...suggestion, author_email: evt.target.value }))
              }}
              leftIcon="envelope"
            />
          </ControlGroup>
          <InputGroup
            type="text"
            placeholder="Affiliation"
            value={suggestion.author_org}
            onChange={evt => {
              setSuggestion(({ ...suggestion }) => ({ ...suggestion, author_org: evt.target.value }))
            }}
            leftIcon="office"
          />
        </ControlGroup>
      </FormGroup>
      <FormGroup
        label="Description"
        labelInfo="(required)"
        helperText={input ? `A description about what this algorithm or data transformation does with the ${input.meta.label || input.spec} with relevant links` : `A description of the core data type`}
      >
        <TextArea
          placeholder={`Your component description goes here`}
          growVertically
          fill
          large
          intent={Intent.PRIMARY}
          onChange={evt => {
            setSuggestion(({ ...suggestion }) => ({ ...suggestion, description: evt.target.value }))
          }}
          value={suggestion.description}
        />
      </FormGroup>
      <Button
        text="Submit"
        disabled={
          !suggestion.name
          || !suggestion.author_name
          || !suggestion.author_email
          || !suggestion.author_org
          || !suggestion.description
          || (!!suggestion.inputs && !suggestion.output)
        }
        onClick={async () => {
          const suggestion_final = {...suggestion}
          if (!suggestion_final.inputs) {
            suggestion_final.name = `Input ${suggestion.name}`
            suggestion_final.output = suggestion.name
          }
          // register the suggestion
          const kvReq = await fetch(`/api/suggest/`, {
            method: 'POST',
            body: JSON.stringify(suggestion_final)
          })
          // construct the KRG node locally
          let OutputNode = krg.getDataNode(suggestion_final.output)
          if (OutputNode === undefined) {
            OutputNode = MetaNode.createData(suggestion_final.output)
              .meta({
                label: suggestion_final.output,
                description: `A data type, suggested as part of ${suggestion_final.name}`,
              })
              .codec<any>()
              .view((props) => {
                return <div>This data type was suggested as part of {suggestion_final.name}</div>
              })
              .build()
            krg.add(OutputNode)
          }
          const ProcessNode = MetaNode.createProcess(suggestion_final.name)
            .meta({
              label: suggestion_final.name,
              description: suggestion_final.description,
            })
            .inputs(dict.init(suggestion_final.inputs.split(',').filter(s => s != '').map((spec, ind) =>
            ({ key: ind.toString(), value: krg.getDataNode(spec) }))))
            .output(OutputNode)
            .prompt((props) => {
              return <div>This was suggested by {suggestion_final.author_name} &lt;{suggestion_final.author_email}&gt; ({suggestion_final.author_org})</div>
            })
            .build()
          krg.add(ProcessNode)
          // extend using those nodes
          const inputs: Record<string, { id: string }> = {}
          if (head) {
            for (const arg in ProcessNode.inputs) {
              inputs[arg] = { id: head.process.id }
            }
          }
          const extendReq = await fetch(`/api/db/fpl/${id}/extend`, {
            method: 'POST',
            body: JSON.stringify({
              type: ProcessNode.spec,
              inputs,
            })
          })
          const extendRes = z.string().parse(await extendReq.json())
          router.push(`/graph/${extendRes}`)
        }}
      />
    </div>
  )
}