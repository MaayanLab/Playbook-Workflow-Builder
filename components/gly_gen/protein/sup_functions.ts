import { URLSearchParams } from "url";
import { z } from "zod";
import {
  GlyGenProteinResponse,
  GlyGenProteinSetResponse,
  ProteinNameAPIResponse,
  GlycosylationAPIResponse,
  GlycosylationData,
  PhosphorylationAPIResponse,
  PhosphorylationData
} from "./data_models";

// Infer the Typescript type from the Zod schemas
type GlyGenProteinResponseType = z.infer<typeof GlyGenProteinResponse>;
type GlyGenProteinSetResponseType = z.infer<typeof GlyGenProteinSetResponse>;
type ProteinNameAPIResponseType = z.infer<typeof ProteinNameAPIResponse>;
type GlycosylationDataType = z.infer<typeof GlycosylationData>;
type GlycosylationAPIResponseType = z.infer<typeof GlycosylationAPIResponse>;
type PhosphorylationDataType = z.infer<typeof PhosphorylationData>;
type PhosphorylationAPIResponseType = z.infer<
  typeof PhosphorylationAPIResponse
>;

export async function get_single_protein_data(
  accession: String,
): Promise<GlyGenProteinResponseType> {
  /* Hits the GlyGen protein detail endpoint and returns the
   * data according to the GlyGenProteinResponse data model.
   *
   * accession: string of the uniprot canonical accession to
   *            search on.
   *
   * Returns: object in the format of GlyGenProteinResponse.
   */

  const detail_response = await fetch(
    `https://api.glygen.org/protein/detail/${accession}/`,
    {
      method: "POST",
      headers: {
        accept: "application/type",
        "Content-Type": "application/json",
      },
    },
  );

  // check if GlyGen API errored (# TODO should probably figure out better error handling to hide
  // from the user)
  if (!detail_response.ok) {
    throw new Error(
      `Protein detail failed with status: ${detail_response.status} ${detail_response.statusText}`,
    );
  }

  const return_data = await detail_response.json();

  // filter on only recommended protein names
  const filtered_protein_names = return_data.protein_names.filter(
    (entry: ProteinNameAPIResponseType) => {
      return entry.type === "recommended";
    },
  );
  return_data.protein_names = filtered_protein_names[0];
  return_data.gene = return_data.gene[0];
  return_data.species = return_data.species[0];
  return_data.species.taxid = return_data.species.taxid.toString();

  // handle glycosylation data
  const is_glycosylation =
    return_data.keywords &&
    return_data.keywords.some((keyword: String) => keyword === "glycoprotein");
  const glycosylation_data = return_data.glycosylation
    ? return_data.glycosylation.map((tag: GlycosylationAPIResponseType) => ({
        site_lbl: tag.site_lbl || "",
        site_category: tag.site_category || "",
        type: tag.type || "",
        glytoucan_ac: tag.glytoucan_ac || "",
      }))
    : [];
  return_data.glycoprotein = {
    glycosylation: is_glycosylation,
    glycosylation_data: glycosylation_data,
    protein_accession: filtered_protein_names[0].name,
  };

  // handle phosphorylation data
  const is_phosphorylation = "phosphorylation" in return_data;
  const phosphorylation_data = return_data.phosphorylation
    ? return_data.phosphorylation.map(
        (tag: PhosphorylationAPIResponseType) => ({
          start_pos: tag.start_pos || "No data available",
          end_pos: tag.end_pos || "No data available",
          kinase_uniprot_canonical_ac:
            tag.kinase_uniprot_canonical_ac || "No data available",
          kinase_gene_name: tag.kinase_gene_name || "No data available",
          residue: tag.residue || "No data available",
          comment: tag.comment || "No data available",
        }),
      )
    : [];
  return_data.phosphorylation = {
    phosphorylation: is_phosphorylation,
    phosphorylation_data: phosphorylation_data,
  };

  const validated_data = GlyGenProteinResponse.parse(return_data);

  return validated_data;
}

export async function get_protein_set_data(accessions: Array<String>) {
  /* Hits the GlyGen protein search endpoint uses the list id
   * to hit the protein list endpoint. If an endpoint is invalid or does
   * not exist, the API fails silently and simply doesn't return data for
   * that accession.
   *
   * accessions: array of uniprot accessions to search on.
   *
   * Returns: array with objects in the format of GlyGenProteinSetResponse.
   */

  // search params for protein search endpoint
  const search_query = new URLSearchParams();
  search_query.append(
    "query",
    JSON.stringify({
      operation: "AND",
      query_type: "search_protein",
      uniprot_canonical_ac: accessions.join(","),
    }),
  );

  const search_query_response = await fetch(
    `https://api.glygen.org/protein/search?${search_query.toString()}`,
    {
      method: "GET",
      headers: {
        accept: "application/json",
        "Content-Type": "application/json",
      },
    },
  );

  // check if GlyGen API errored (# TODO should probably figure out better error handling to hide
  // from the user)
  if (!search_query_response.ok) {
    throw new Error(
      `Protein search endpoint failed with status: ${search_query_response.status} ${search_query_response.statusText}`,
    );
  }
  const search_query_result = await search_query_response.json();
  const list_id = search_query_result.list_id;

  // search params for protein list endpoint
  const list_query = new URLSearchParams();
  list_query.append(
    "query",
    JSON.stringify({
      id: list_id,
      offset: 1,
      limit: accessions.length,
      order: "desc",
      sort: "hit_score",
      filters: [],
    }),
  );
  const list_query_response = await fetch(
    `https://api.glygen.org/protein/list?${list_query.toString()}`,
    {
      method: "GET",
      headers: {
        accept: "application/json",
        "Content-Type": "application/json",
      },
    },
  );

  // check if GlyGen API errored (# TODO should probably figure out better error handling to hide
  // from the user)
  if (!list_query_response.ok) {
    throw new Error(
      `Protein search endpoint failed with status: ${search_query_response.status} ${search_query_response.statusText}`,
    );
  }
  const list_query_result = await list_query_response.json();
  const list_query_results = list_query_result.results;
  const return_data = [];
  for (const element of list_query_results) {
    return_data.push({
      gene: {
        name: element.gene_name,
      },
      uniprot: {
        uniprot_canonical_ac: element.uniprot_canonical_ac,
      },
      protein_names: {
        name: element.protein_name,
      },
      species: {
        name: element.organism,
        taxid: String(element.tax_id),
      },
      bools: {
        total_n_glycosites: element.total_n_glycosites,
        total_o_glycosites: element.total_o_glycosites,
        reported_phosphosites: element.reported_phosphosites,
        reported_snv: element.reported_snv,
      },
    });
  }

  return return_data;
}

export function extract_specific_data(data: GlyGenProteinResponseType, key: String): GlycosylationDataType | PhosphorylationDataType {
  /* Takes the full protein response and extracts specific data such as the 
   * glycosylation or phosphorylation data. 
   *
   * data: the full protein detail API response. 
   * key: the data key to extract (accepts glycosylation or phosphorylation).
   *
   * Returns: the extracted object. 
   */
  
  if (key.toLowerCase().trim() === "glycosylation") {
    return data.glycoprotein;
  } else if (key.toLowerCase().trim() === "phosphorylation") {
    return data.phosphorylation;
  } else {
    throw new Error("Invalid key type to extract_specific_data.")
  }
}

// --- Type Guards ---//

export function glycosylation_check(data: any): data is GlycosylationDataType {
  return (data as GlycosylationDataType).glycosylation !== undefined;
}

export function phosphorylation_check(data: any): data is PhosphorylationDataType {
  return (data as PhosphorylationDataType).phosphorylation !== undefined;
}
