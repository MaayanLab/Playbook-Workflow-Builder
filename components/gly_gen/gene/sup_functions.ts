import { get_single_protein_data } from "../protein/sup_functions";
import { GlyGenProteinResponse, ProteinNameAPIResponse } from "../protein/data_models";
import { z } from "zod";

// Infer the Typescript type from Zod schemas 
type GlyGenProteinResponseType = z.infer<typeof GlyGenProteinResponse>;

export async function filter_glygen_proteins(gene_symbol: String): Promise<GlyGenProteinResponseType> {
  /* Takes a gene symbol and finds the protein products from GlyGen.
   *
   * gene_symbol: The gene symbol to search on.
   *
   * Returns: object in the format of GlyGenProteinResponse.
   */

  const search_response = await fetch("https://api.glyGen.org/protein/search/", {
    method: "POST",
    headers: {
      accept: "application/json",
      "Content-Type": "application/json",
    },
    body: JSON.stringify({gene_name: gene_symbol})
  })

  // check if GlyGen API errored (# TODO should probably figure out better error handling to hide
  // from the user)
  if (!search_response.ok) {
    throw new Error(`Protein detail failed with status: ${search_response.status} ${search_response.statusText}`,)
  };
  
  const search_result = await search_response.json();
  const id = search_result.list_id;
  
  const protein_response = await fetch(
    "https://api.glyGen.org/protein/list/",
    {
      method: "POST",
      headers: {
        accept: "application/json",
        "Content-Type": "application/json"
      },
      body: JSON.stringify({id: id})
    }
  );

  // check if GlyGen API errored (# TODO should probably figure out better error handling to hide
  // from the user)
  if (!protein_response.ok) {
    throw new Error(`Protein detail failed with status: ${protein_response.status} ${protein_response.statusText}`,)
  };

  const protein_result = await protein_response.json();

  // filter list results 
  console.log("Gene symbol: %s", gene_symbol)
  for (const item of protein_result["results"]) {
    if (item.hasOwnProperty("gene_name")) {
      const gene_name = item["gene_name"].toUpperCase();
      const species = item["organism"];
      console.log("----------------------------")
      console.log("-> Prop name: %s", gene_name);
      console.log("-> Species: %s", species);
      if (gene_name === gene_symbol.toUpperCase() && species === "Human") {
        console.log("==========================")
        console.log(`Human result: ${item["gene_name"]}`)
        console.log(`Type of result: ${typeof(item)}`)
        console.log("==========================")
        return get_single_protein_data(item["uniprot_canonical_ac"]);
      }
    }
  }
  throw new Error("No matching protein product found.");
}
