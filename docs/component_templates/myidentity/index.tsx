import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/term'

// A unique name for your resolver is used here
export const MyIdentity = MetaNode('MyIdentity')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'My Identity Function',
    description: 'Compute the identity function',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ input: GeneTerm })
  // This should be a single Meta Node Data Type
  .output(GeneTerm)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    return props.inputs.input
  })
  .story(props => ({ abstract: `The identity function is applied to${props.inputs ? ` ${props.inputs.input}` : ''}.` }))
  .build()
