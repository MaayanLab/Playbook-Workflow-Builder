export function filterGlyGenResults(result, gene_name) {
    for (const idx in result['results']) {
      if (result['results'][idx].hasOwnProperty('gene_name')){
        const geneName = result['results'][idx]['gene_name'].toUpperCase();
        const speciesName = result['results'][idx]['organism'];
        if (geneName === gene_name.toUpperCase() && speciesName === 'Homo sapiens'){ 
            console.log('==========================');
            console.log(`Human result: ${result['results'][idx]['gene_name']}`)
            console.log('==========================');
            return result['results'][idx]
        }
      }
    }
  }

