import React from 'react'
import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/term'

// A unique name for your prompt is used here
export const MyPrompt = MetaNode('MyPrompt')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'My Prompt',
    description: 'Prompt the user',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  // Prompts don't always call for inputs, thus it can be blank.
  .inputs()
  // This must be a single Meta Node Data Type
  .output(GeneTerm)
  // The prompt function is a react component which the user will
  //  engage with. The `submit` function should be used to resolve
  //  the prompt, ultimately providing data in the shape prescribed
  //  by the data type codec of the output.
  .prompt(props => {
    const [gene, setGene] = React.useState('')
    // the output parameter of a prompt will be present if the
    //  prompt has previously been populated, and it may modify
    //  without remount, using useEffect to update state with the
    //  output when it changes is ideal.
    React.useEffect(() => {
      setGene(props.output || '')
    }, [props.output])
    return (
      <div>
        <input value={gene} onChange={evt => setGene(evt.target.value)} />
        <button onClick={evt => {
          // when the user interaction has resolved, we use the
          //  submit callback to resolve this prompt
          props.submit(gene)
        }}>Submit</button>
      </div>
    )
  })
  // Write up what this node has done doing
  .story(props => `The workflow starts with ${props.output ? props.output : 'a gene'}.`)
  .build()
