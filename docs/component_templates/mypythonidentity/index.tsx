import python from '@/utils/python'
import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/input/term'

// A unique name for your resolver is used here
export const MyPythonIdentity = MetaNode('MyPythonIdentity')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'My Python Identity Function',
    description: 'Compute the identity function',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ input: GeneTerm })
  // This should be a single Meta Node Data Type
  .output(GeneTerm)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => await python(
    // this is the absolute import location of the function
    //  we'll call during the resolution phase. it is important
    //  to make sure the python function returns a datatype which
    //  is json serializable in a form to match the output data type codec
    'components.mypythonidentity.identity',
    // these kargs/kwargs are passed to the python function
    //  i.e. identity(*kargs, **kwargs)
    { kargs: [props.inputs.input], kwargs: {} },
  ))
  .build()
