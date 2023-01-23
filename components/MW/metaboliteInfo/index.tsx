import { MetaNode } from '@/spec/metanode'
import { MetaboliteTerm } from '@/components/core/input/term'
import { MetaboliteSet } from '@/components/core/input/set'
import { MetaboliteSummary } from '../metabolite_summary'
// A unique name for your resolver is used here
export const MetaboliteInfo = MetaNode.createProcess('MetaboliteInfo')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Metabolite information',
    description: 'Extract information for a metabolite',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ metabolite: MetaboliteTerm })
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

  export const MetaboliteSetInfo = MetaNode.createProcess('MetaboliteSetInfo')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Metabolite set information',
    description: 'Extract information for a set of metabolites',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ metaboliteset: MetaboliteSet })
  // This should be a single Meta Node Data Type
  .output(MetaboliteSummary)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    const metNameSet = props.inputs.metaboliteset
    var jsonArrayObject = []
    //metNameSet.forEach(async metName => {
    for (let metName of metNameSet) {
      const req = await fetch(`https://www.metabolomicsworkbench.org/rest/refmet/name/${metName}/all`)
      const res = await req.json()
      //console.log(res)
      jsonArrayObject.push(res)
    }
    
    await console.log(jsonArrayObject)
    
    return jsonArrayObject
    //return props.inputs.input
  })
  .build()