import { MetaNode } from '@/spec/metanode'
import { GeneTerm } from '@/components/core/term'
import { GeneSet } from '@/components/core/set'
import { ConvertedGeneID } from '@/components/MW/ConvertedGeneID'
import { uniqJsonSubset } from '@/components/MW/ConvertedGeneID'
import { metabolomicsworkbench_icon } from '@/icons'

// A unique name for your resolver is used here
export const GeneIDConv = MetaNode('GeneIDConv')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Gene ID Conversion',
    description: 'Given one type of gene ID, generate other types of gene IDs',
    icon: [metabolomicsworkbench_icon],
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ gene: GeneTerm })
  // This should be a single Meta Node Data Type
  .output(ConvertedGeneID)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    // example 1: https://bdcw.org/geneid/rest/species/hsa/GeneIDType/SYMBOL_OR_ALIAS/GeneListStr/ITPR3__IL6__KLF4/View/json
    // example 2: https://bdcw.org/geneid/rest/species/hsa/GeneIDType/SYMBOL/GeneListStr/ITPR3__IL6__KLF4/View/json
    // example 3: https://bdcw.org/geneid/rest/species/hsa/GeneIDType/ENTREZID/GeneListStr/3569,3710,9314/View/json
    const species_id = "hsa"
    const geneid_type = "SYMBOL_OR_ALIAS"
    const gene_id = props.inputs.gene // "HK1" //"${props.inputs.gene}"
    const json_or_txt = "json"
    const req = await fetch(`https://bdcw.org/geneid/rest/species/${species_id}/GeneIDType/${geneid_type}/GeneListStr/${gene_id}/View/${json_or_txt}`)
    const res = await req.json()

    //return props.inputs.gene
    return res
  })
  .story(props => ({
    abstract: `The gene ID was then convert to various other alternative gene IDs using the Metabolomics Gene Conversion API\\ref{Gene ID Conversion Tool, https://bdcw.org/geneid/}.`
  }))
  .build()

// For GeneSet
export const GeneSetIDConv = MetaNode('GeneSetIDConv')
  // Human readble descriptors about this node should go here
  .meta({
    label: 'Gene ID Conversion for a set of genes',
    description: 'Given one type of gene ID for a set of genes, generate other types of gene IDs.',
  })
  // This should be a mapping from argument name to argument type
  //  the types are previously defined Meta Node Data Types
  .inputs({ geneset : GeneSet  })
  // This should be a single Meta Node Data Type
  .output(ConvertedGeneID)
  // The resolve function uses the inputs and returns output
  //  both in the shape prescribed by the data type codecs
  .resolve(async (props) => {
    // example 1: https://bdcw.org/geneid/rest/species/hsa/GeneIDType/SYMBOL_OR_ALIAS/GeneListStr/ITPR3__IL6__KLF4/View/json
    // example 2: https://bdcw.org/geneid/rest/species/hsa/GeneIDType/SYMBOL/GeneListStr/ITPR3__IL6__KLF4/View/json
    // example 3: https://bdcw.org/geneid/rest/species/hsa/GeneIDType/ENTREZID/GeneListStr/3569,3710,9314/View/json
    const species_id = "hsa"
    const geneid_type = "SYMBOL_OR_ALIAS"
    const gene_id = props.inputs.geneset.set.join(","); // "HK1" //"${props.inputs.gene}"
    //console.log(gene_id)
    const json_or_txt = "json"
    const req = await fetch(`https://bdcw.org/geneid/rest/species/${species_id}/GeneIDType/${geneid_type}/GeneListStr/${gene_id}/View/${json_or_txt}`)
    const res = await req.json()

    //return props.inputs.gene
    return res
  })
  .story(props => ({
    abstract: `The gene IDs were converted to various other alternative gene IDs using the Metabolomics Gene Conversion API\\ref{Gene ID Conversion Tool, https://bdcw.org/geneid/}.`
  }))
  .build()

// Process to convert ConvertedGeneID to GeneInfo
// Much of LINCS's APIs are centered around GeneInfo
export const ConvertedGeneID2GeneInfo = MetaNode('ConvertedGeneID2GeneInfo')
  .meta({
    label: 'Extract Gene Term from ConvertedGeneID',
    description: 'Given a ConvertedGeneID object, convert it to Gene Term object',
  })
  .inputs({data: ConvertedGeneID})
  .output(GeneTerm)
  .resolve(async (props) => {
    return uniqJsonSubset(props.inputs.data)[0].SYMBOL;
  })
  .story(props => ({
    abstract: `The official gene symbol${props.output ? ` ${props.output}` : ''} was resolved.`
  }))
  .build()
