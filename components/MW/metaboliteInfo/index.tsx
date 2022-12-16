import { MetaNode } from '@/spec/metanode'
import { MetaboliteName } from '@/components/MW/metabolite'
import { MetaboliteSummary } from '../metabolite_summary'
// A unique name for your resolver is used here
export const MetaboliteInfo = MetaNode.createProcess('MetaboliteInfo')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'My Metabolite Info Function',
    description: 'Compute the metabolite info function',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ metabolite: MetaboliteName })
  // This should be a single Meta Node Data Type
  .output(MetaboliteSummary)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const metName = props.inputs.metabolite
    const req = await fetch(`https://www.metabolomicsworkbench.org/rest/refmet/name/${metName}/all`)
    const res = await req.json()
    
    return  res
    //return props.inputs.input
  })
  .build()
