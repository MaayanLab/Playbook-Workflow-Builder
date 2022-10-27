import React from 'react'
import { MetaNode, MetaNodeDataType } from '@/spec/metanode'
import { z } from 'zod'
import { Intent } from '@blueprintjs/core'
import dynamic from 'next/dynamic'

const Button = dynamic(() => import('@blueprintjs/core').then(({ Button }) => Button))
const FormGroup = dynamic(() => import('@blueprintjs/core').then(({ FormGroup }) => FormGroup))
const ControlGroup = dynamic(() => import('@blueprintjs/core').then(({ ControlGroup }) => ControlGroup))
const InputGroup = dynamic(() => import('@blueprintjs/core').then(({ InputGroup }) => InputGroup))
const TextArea = dynamic(() => import('@blueprintjs/core').then(({ TextArea }) => TextArea))
const MenuItem = dynamic(() => import('@blueprintjs/core').then(({ MenuItem }) => MenuItem))
const Suggest2 = dynamic(() => import('@blueprintjs/select').then(({ Suggest2 }) => Suggest2))

export const Suggestion = MetaNode.createData('Suggestion')
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

export const SuggestInputEdge = MetaNode.createProcess('SuggestInputEdge')
  .meta({
    label: 'Suggest an input type',
    description: 'Provide a description about what should be here',
  })
  .inputs()
  .output(Suggestion)
  .prompt((props) => {
    const [suggestion, setSuggestion] = React.useState('')
    return (
      <TextArea
        growVertically={true}
        large={true}
        intent={Intent.PRIMARY}
        onChange={(evt) => setSuggestion(evt.target.value)}
        value={suggestion}
      />
    )
  })
  .build()


import { FileURL } from '@/components/core/file'
import { GeneCountMatrix } from '@/components/data/gene_count_matrix'
import { GeneTerm } from '@/components/core/input/term'
import { GeneInfo } from '@/components/service/mygeneinfo'
import { SignificantTissues } from '@/components/core/significant_tissues'
import krg from '@/app/krg'

export const SuggestInteractiveEdge = [
  FileURL,
  GeneCountMatrix,
  GeneTerm,
  GeneInfo,
  SignificantTissues,
].flatMap(input => [
  MetaNode.createProcess(`SuggestInteractiveEdge[${input.spec}]`)
    .meta({
      label: 'Suggest a visualization method',
      description: `This would visualize the ${input.meta.label || input.spec}. Provide a description about what should be here.`,
    })
    .inputs({ input } as Record<string, MetaNodeDataType>)
    .output(Suggestion)
    .prompt((props) => {
      const [suggestion, setSuggestion] = React.useState('')
      return (
        <TextArea
          growVertically={true}
          large={true}
          intent={Intent.PRIMARY}
          onChange={(evt) => setSuggestion(evt.target.value)}
          value={suggestion}
        />
      )
    })
    .build(),
  MetaNode.createProcess(`SuggestResolveEdge[${input.spec}]`)
    .meta({
      label: 'Suggest an algorithm or data transformation method',
      description: `This would transform the ${input.meta.label || input.spec}. Provide a description about what should be here.`,
    })
    .inputs({ input } as Record<string, MetaNodeDataType>)
    .output(Suggestion)
    .prompt((props) => {
      const [suggestion, setSuggestion] = React.useState({
        name: '',
        inputs: input.spec as string,
        output: '',
        author_name: '',
        author_email: '',
        author_org: '',
        description: '',
      })
      React.useEffect(() => {
        if (props.output) setSuggestion(props.output)
      }, [props.output])
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
            helperText={`A description about what this algorithm or data transformation does with the ${input.meta.label || input.spec} with relevant links`}
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
            onClick={() => props.submit(suggestion)}
          />
        </div>
      )
    })
    .build()
  ]
)
