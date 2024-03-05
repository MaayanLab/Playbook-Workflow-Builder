import { URLSearchParams } from "url";
import { z } from "zod";

export async function glycan_set_search_query(glytoucan_accessions: String[]) {
  /* Retrieves the Glycan set data from the GlyGen APIs.
   *
   * glytoucan_accessions: Array glycan glytoucan accessions.
   *
   * Returns: Glycan set object.
   */
  const search_query = new URLSearchParams();
  search_query.append(
    "query",
    JSON.stringify({
      operation: "AND",
      query_type: "search_glycan",
      glycan_identifier: {
        glycan_id: glytoucan_accessions.join(","),
      },
    }),
  );

  const search_query_response = await fetch(
    `https://api.glygen.org/glycan/search?${search_query.toString()}`,
    {
      method: "GET",
      headers: {
        accept: "application/json",
        "Content-Type": "aapplication/json",
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

  const list_query = new URLSearchParams();
  list_query.append(
    "query",
    JSON.stringify({
      id: list_id,
      offset: 1,
      limit: glytoucan_accessions.length,
      order: "desc",
      sort: "hit_score",
      filters: [],
    }),
  );
  const list_query_response = await fetch(
    `https://api.glygen.org/glycan/list?${list_query.toString()}`,
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
      `Protein search endpoint failed with status: ${list_query_response.status} ${list_query_response.statusText}`,
    );
  }
  const list_query_result = await list_query_response.json();
  const list_query_results = list_query_result.results;
  const result = [];
  for (const element of list_query_results) {
    result.push({
      glytoucan: {
        glytoucan_ac: element.glytoucan_ac,
      },
      hit_score: element.hit_score,
      mass: element.mass,
      mass_pme: element.mass_pme,
      sugar_count: element.number_monosaccharides,
      glycoprotein_count: element.number_proteins,
      associated_enzymes: element.number_enzymes,
    });
  }

  return result;
}
