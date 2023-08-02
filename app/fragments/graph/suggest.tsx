import { NextRouter, useRouter } from 'next/router'
import type { Metapath } from '@/app/fragments/metapath'
import React from 'react'
import { MetaNode, DataMetaNode, ProcessMetaNode } from '@/spec/metanode'
import { z } from 'zod'
import dynamic from 'next/dynamic'
import * as dict from '@/utils/dict'
import type KRG from '@/core/KRG'
import { z_uuid } from '@/utils/zod'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'

const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))
const FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const ControlGroup = dynamic(() => import('@blueprintjs/core').then(({ ControlGroup }) => ControlGroup))
const InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))
const TextArea = dynamic(() => import('@blueprintjs/core').then(({ TextArea }) => TextArea))
const MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))
const Suggest2 = dynamic(() => import('@blueprintjs/select').then(({ Suggest2 }) => Suggest2))
const UserIdentity = dynamic(() => import('@/app/fragments/graph/useridentity'))

const Suggestion = MetaNode('Suggestion')
  .meta({
    label: 'Suggestion',
    description: 'An actual suggestion',
  })
  .codec(z.object({
    name: z.string(),
    inputs: z.string(),
    output: z.string(),
    user: z_uuid(),
    description: z.string()
  }))
  .view(suggestion => (
    <div>
      <div className="bp4-card">
        <h4 className="bp4-heading">{suggestion.name}</h4>
        <p className="bp4-text-large">{suggestion.description}</p>
        <div className="bp4-callout">
          <h5 className="bp4-heading">Author</h5>
          <UserIdentity user={suggestion.user} />
        </div>
        {/* <br />
        <h5 className="bp4-heading"><i>Comments</i></h5>
        <div className="bp4-callout">
          <p className="bp4-text-large">My comment is important</p>
          <h5 className="bp4-heading">Author</h5>
          <UserIdentity user={suggestion.user} />
        </div> */}
      </div>
    </div>
  ))
  .build()

export function SuggestionEdges(input?: DataMetaNode) {
  const suggestion_edges: Array<ProcessMetaNode> = []
  if (input === undefined) {
    suggestion_edges.push(
      MetaNode(`SuggestDataType`)
        .meta({
          label: 'Suggest a core data type',
          description: `This would be usable as an initial or intermediary data`,
          pagerank: -99,
        })
        .inputs()
        .output(Suggestion)
        .prompt(props => <></>)
        .story(props => ``)
        .build() as ProcessMetaNode
    )
  }
  if (input !== undefined) {
    suggestion_edges.push(
      MetaNode(`SuggestInteractiveEdge[${input.spec}]`)
        .meta({
          label: 'Suggest a visualization method',
          description: `This would visualize the ${input.meta.label || input.spec}. Provide a description about what should be here.`,
          pagerank: -99,
        })
        .inputs({ input })
        .output(Suggestion)
        .prompt((props) => <></>)
        .story(props => ``)
        .build() as ProcessMetaNode
    )
    suggestion_edges.push(
      MetaNode(`SuggestResolveEdge[${input.spec}]`)
        .meta({
          label: 'Suggest an algorithm or data transformation method',
          description: `This would transform the ${input.meta.label || input.spec}. Provide a description about what should be here.`,
          pagerank: -99,
        })
        .inputs({ input })
        .output(Suggestion)
        .prompt((props) => <></>)
        .story(props => ``)
        .build() as ProcessMetaNode
    )
  }
  return suggestion_edges.map(element => ({
    ...element,
    onClick: ({ router, id, head }: { router: NextRouter, id: string, head: Metapath }) => {
      router.push(`/graph/${id}${head ? head.id !== id ? `/node/${head.id}` : '' : '/node/start'}/suggest`)
    }
  }))
}

export default function Suggest({ session_id, krg, id, head }: { session_id?: string, krg: KRG, id: string, head: Metapath }) {
  const router = useRouter()
  const { data: userSession } = useSessionWithId({ required: true })
  const processNode = head ? krg.getProcessNode(head.process.type) : undefined
  const input = processNode ? processNode.output : undefined
  const [suggestion, setSuggestion] = React.useState({
    name: '',
    inputs: input ? input.spec as string : '',
    output: '',
    user: userSession?.user?.id,
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
        label="Description"
        labelInfo="(required)"
        helperText={input ? `A description about what this algorithm or data transformation does with the ${input.meta.label || input.spec} with relevant links` : `A description of the core data type`}
      >
        <TextArea
          placeholder={`Your component description goes here`}
          growVertically
          fill
          large
          intent="primary"
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
          || !suggestion.user
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
            OutputNode = MetaNode(suggestion_final.output)
              .meta({
                label: `${suggestion_final.output} (Suggestion)`,
                description: `A data type, suggested as part of ${suggestion_final.name}`,
                pagerank: -100,
              })
              .codec<any>()
              .view((props) => {
                return <div>This data type was suggested as part of {suggestion_final.name}</div>
              })
              .build()
            krg.add(OutputNode)
          }
          const ProcessNode = MetaNode(suggestion_final.name)
            .meta({
              label: `${suggestion_final.name} (Suggestion)`,
              description: suggestion_final.description,
              pagerank: -100,
            })
            .inputs(dict.init(suggestion_final.inputs.split(',').filter(s => s != '').map((spec, ind) =>
            ({ key: ind.toString(), value: krg.getDataNode(spec) })).filter(({ key, value }) => !!value)))
            .output(OutputNode)
            .prompt((props) => {
              return <div>
                <p>{suggestion.description}</p>
                <p>This was suggested by {suggestion.user ? <UserIdentity user={suggestion.user} /> : <>a playbook partnership user</>}.</p>
              </div>
            })
            .story(props => `It is suggested that "${suggestion.description}" be applied to the inputs: ${suggestion.inputs} to get a ${OutputNode.meta.label}.`)
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