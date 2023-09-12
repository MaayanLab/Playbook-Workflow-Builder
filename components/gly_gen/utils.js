// import * as fs from 'fs';

async function resolveFilteredResult(cannonicalAccession) {
  const resolvedDetail = await fetch(`https://api.glygen.org/protein/detail/${cannonicalAccession}/`, {
    method: 'POST',
    headers: {
      accept: 'application/json',
     'Content-Type': 'application/json',
    },
  });
  const intermediateResult = await resolvedDetail.json();
  // fs.writeFileSync('intermediate_result.json', JSON.stringify(intermediateResult, null, 2));
  const filteredProteinName = intermediateResult.protein_names.filter(
    (name) => { 
      return(name.type === "recommended"); 
    }
  );
  // TODO: Is there ever a case where gene-as-array is appropriate?
  intermediateResult.gene = intermediateResult.gene[0];
  // TODO: Is there ever a case where species-as-array is appropriate?
  intermediateResult.species = intermediateResult.species[0];
  intermediateResult.protein_names = filteredProteinName[0];
  intermediateResult.species.taxid = intermediateResult.species.taxid.toString();
  // Q96F25: false test case 
  // HGF: true test case
  intermediateResult.glycoprotein = {
    glycosylation: intermediateResult.keywords && intermediateResult.keywords[0] === 'glycoprotein'
  }
  return intermediateResult;
}

export function filterGlyGenResults(result, prop_type, prop_name) {
    
    const result_type = prop_type === 'gene' ? 'gene_name' : 'protein_name'

    for (const idx in result['results']) {
      if (result['results'][idx].hasOwnProperty(result_type)){
        const propName = result['results'][idx][result_type].toUpperCase();
        const speciesName = result['results'][idx]['organism'];
        console.log('-> Prop name: %s', propName);
        console.log('-> Species name: %s', speciesName);
        console.log('+> Species match?: %s', propName === prop_name.toUpperCase() && speciesName === 'Homo sapiens')
        if (propName === prop_name.toUpperCase() && speciesName === 'Homo sapiens'){
            console.log('==========================');
            console.log(`Human result: ${result['results'][idx]['gene_name']}`)
            console.log(`Type of result: ${typeof(result['results'][idx])}`)
            console.log('==========================');
            return resolveFilteredResult(result['results'][idx]['uniprot_canonical_ac']);
        }
      }
    }
  }

