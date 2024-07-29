import { MetaNode } from '@/spec/metanode'
import { MetaboliteTerm } from '@/components/core/term'
import { MetaboliteSet } from '@/components/core/set'
import { MetaboliteSummary } from '@/components/MW/metabolite_summary'
import { metabolomicsworkbench_icon } from '@/icons'

// A unique name for your resolver is used here
export const MetaboliteInfo = MetaNode('MetaboliteInfo')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Metabolite information',
    description: 'Extract information for a metabolite',
    icon: [metabolomicsworkbench_icon],
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ metabolite: MetaboliteTerm })
  // This should be a single Meta Node Data Type
  .output(MetaboliteSummary)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    var jsonArrayObject = []; // Added 2023/01/25
    const metName = props.inputs.metabolite
    //Original by Sumana: const req = await fetch(`https://www.metabolomicsworkbench.org/rest/refmet/name/${metName}/all`)
    // Mano: 2023/08/01: using 'match' instead of 'name' so that it can take even kegg_id
    const req = await fetch(`https://www.metabolomicsworkbench.org/rest/refmet/match/${metName}/all`)    
    const res = await req.json()

    //return  res ; // older

    jsonArrayObject.push(res); // Added 2023/01/25

    return jsonArrayObject; // Added 2023/01/25
    //return props.inputs.input
  })
  .story(props => ({
    abstract: `The metabolite was then searched in the Metabolomics Workbench\\ref{The Metabolomics Workbench, https://www.metabolomicsworkbench.org/} to extract more information about the metabolite.`,
    introduction: `The Metabolomics Workbench Metabolite Database contains structures and annotations of biologically relevant metabolites. MW supports REST APIs that provides information regarding the metabolite given the metabolite name.  `,
    methods: `Given a metabolite name, say, metName, the REST API in MW,https://www.metabolomicsworkbench.org/rest/refmet/match/${props.inputs?.metabolite ? props.inputs.metabolite : '{metName}'}/all, provides the metabolite information in JSON format that is displayed in HTML.`,
    legend: `The metabolite information summarizes the refmet name, formula, exact mass, super, main and sub classes of the metabolite.`,
  }))
  .build()

export const MetaboliteSetInfo = MetaNode('MetaboliteSetInfo')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Extract metabolite set information',
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
    const metNameSet = props.inputs.metaboliteset.set
    var jsonArrayObject = []
    //metNameSet.forEach(async metName => {
    for (let metName of metNameSet) {
      //Original by Sumana: const req = await fetch(`https://www.metabolomicsworkbench.org/rest/refmet/name/${metName}/all`)
      // Mano: 2023/08/01: using 'match' instead of 'name' so that it can take even kegg_id
      const req = await fetch(`https://www.metabolomicsworkbench.org/rest/refmet/match/${metName}/all`)

    const res = await req.json()

      //console.log(res)
      // Mano: 2023/06/28: if res is null or [] (for metabolites like h+, CO2, H2O), it gets pushed as element and 
      // creates issue in viewing/rendering MetaboliteSummary in the file ../metabolite_summary/index.tsx
      // If so, may be just create dummy object with only metName
      const res_dummy = {        name: metName,        pubchem_cid: "",        inchi_key: "",
        exactmass: "",        formula: "",        super_class: "",        main_class: "",        sub_class: "",
      }
      if(!Array.isArray(res)){ // if(Array.isArray(res) && res.length == 0)
        jsonArrayObject.push(res)
      } else{
        jsonArrayObject.push(res_dummy)
      }
    }

    //await console.log(jsonArrayObject)

    return jsonArrayObject
    //return props.inputs.input
  })
  .story(props => ({
    abstract: `The metabolites were then searched in the Metabolomics Workbench\\ref{The Metabolomics Workbench, https://www.metabolomicsworkbench.org/} to extarct more information about the metabolites.`,
    introduction: `The Metabolomics Workbench Metabolite Database contains structures and annotations of biologically relevant metabolites. MW supports REST APIs that provides information regarding the metabolite given the metabolite name.`,
   methods: `For each metabolite in the meatbolite set, the metabolite information is obtained from the MW REST API and the resultant JSON is displayed in a tabular format.`,
   legend: `For each metabolite in the set of metabolites, the metabolite information metanode, summarizes the refmet name, formula, exact mass, and the names of super , main and subclasses , the metabolite belongs.`,
  }))
  .build()