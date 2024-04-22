import React from "react";
import { MetaNode } from "@/spec/metanode";
import {
  GeneInfo,
  GeneInfoFromGeneTerm,
} from "@/components/service/mygeneinfo";
import { z } from "zod";
import { glygen_icon } from "@/icons";
import {
  GeneTerm,
  ProteinTerm,
  GlycanTerm,
} from "@/components/core/term";
import { GlyGenProteinResponseNode } from "../protein";
import { filter_glygen_proteins } from "./sup_functions";

// --- Process Metanodes --- //

// Protein process metanode for searching the GlyGen database by gene name for protein
// products given a GeneInfo
export const GlyGenProteinProduct = MetaNode("GGPP")
  .meta({
    label: "Search GlyGen by Gene Name for Protein Products",
    description: "Find protein product records in GlyGen for the gene",
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ gene: GeneInfo })
  .output(GlyGenProteinResponseNode)
  .resolve(async (props) => {
    return filter_glygen_proteins(props.inputs.gene.symbol);
  })
  .story(
    (props) =>
      `Next, the GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of proteins that originate from ${props.inputs ? props.inputs.gene.symbol : "the gene"}.`,
  )
  .build();

// Protein process metanode that resolves a GeneTerm to a protein
export const GlyGenProteinInformation = MetaNode("GlyGenProteinInformation")
  .meta({
    label: "Search GlyGen for Protein Products",
    description: "Find protein product records in GlyGen for the gene",
    icon: [glygen_icon],
    pagerank: 2,
  })
  .inputs({ gene: GeneTerm })
  .output(GlyGenProteinResponseNode)
  .resolve(async (props) => {
    const gene = await GeneInfoFromGeneTerm.resolve(props);
    return await GlyGenProteinProduct.resolve({ ...props, inputs: { gene } });
  })
  .story(
    (props) =>
      `The GlyGen database [\\ref{doi:10.1093/glycob/cwz080}] was searched to identify a relevant set of protein products that originate from ${props.inputs ? props.inputs.gene : "the gene"}.`,
  )
  .build();
