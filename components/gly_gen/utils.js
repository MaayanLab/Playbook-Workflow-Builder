async function resolveFilteredResult(cannonicalAccession) {
  const resolvedDetail = await fetch(`https://api.glygen.org/protein/detail/${cannonicalAccession}/`, {
    method: 'POST',
    headers: {
      accept: 'application/json',
     'Content-Type': 'application/json',
    },
  });
  const intermediateResult = await resolvedDetail.json();
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
  return intermediateResult;
}

export function filterGlyGenResults(result, gene_name) {
    for (const idx in result['results']) {
      if (result['results'][idx].hasOwnProperty('gene_name')){
        const geneName = result['results'][idx]['gene_name'].toUpperCase();
        const speciesName = result['results'][idx]['organism'];
        if (geneName === gene_name.toUpperCase() && speciesName === 'Homo sapiens'){ 
            console.log('==========================');
            console.log(`Human result: ${result['results'][idx]['gene_name']}`)
            console.log(`Type of result: ${typeof(result['results'][idx])}`)
            console.log('==========================');
            return resolveFilteredResult(result['results'][idx]['uniprot_canonical_ac']);
        }
      }
    }
  }

