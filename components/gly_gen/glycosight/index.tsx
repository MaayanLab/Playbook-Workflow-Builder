import React from "react";
import { MetaNode } from "@/spec/metanode";
import { z } from "zod";
import { Table, Cell, Column, EditableCell } from '@/app/components/Table';
import * as dict from '@/utils/dict'
import { downloadBlob } from '@/utils/download'

import { file_icon, file_transfer_icon, glygen_icon } from "@/icons";
import { FileC, FilePrompt, FileURL, FileInput } from "@/components/core/file";
import { fileAsStream } from "@/components/core/file/api/download";
import SafeRender from "@/utils/saferender";

import { GlycoSightOutputC, 
         GlycoSightFileURLC, 
         GlycoSightOutputType,
         GlycoSightFileURLType } from "./data_models";
import { TestGlycoSightAPI,
         UploadAndAnalyze
 } from "./sup_functions";
import { clientLoadExample } from "@/components/gly_gen/api/glycosight/client";
import { notify_insertion_trigger } from "@/db";

export const GlycoSightOutputNode = MetaNode("GlycoSightOutputNode")
    .meta({
        label: `GlycoSight Analysis output`,
        description: `UniProt accessions, peptides, and associated glycosylation sites`,
        icon: [glygen_icon],
    })
    .codec(GlycoSightOutputC)
    .view((data) => {
        const gs_items = data.results;
        // const gs_items = dict.items(data);
        return (
              <Table
                height={500}
                cellRendererDependencies={[gs_items]}
                numRows={gs_items.length}
                enableGhostCells
                enableFocusedCell
                downloads={{
                  JSON: () => downloadBlob(new Blob([JSON.stringify(gs_items)], { type: 'application/json;charset=utf-8' }), 'data.json'),
                //   GS_TSV: () => downloadBlob(new Blob([gs_items.map(({ key: term, value: { description, set } }) => [term, description||'', ...set].join('\t')).join('\n')], { type: 'text/tab-separated-values;charset=utf-8' }), 'data.gs'),
                }}
              >
                {/* <Column
                  name="Index"
                  cellRenderer={row => <Cell key={row+''}>{gs_items[row].key.toString()}</Cell>}
                /> */}
                <Column
                  name="UniProt Accession"
                  cellRenderer={row => <Cell key={row+''}>{gs_items[row].UniProtAcc}</Cell>}
                />
                <Column
                  name="Amino Acid Position"
                  cellRenderer={row => <Cell key={row+''}>{gs_items[row].AAPosition}</Cell>}
                />
                <Column
                  name="Spectral Count"
                  cellRenderer={row => <Cell key={row+''}>{gs_items[row].SpectralCount}</Cell>}
                />
                <Column
                  name="Distinct Peptide Count"
                  cellRenderer={row => <Cell key={row+''}>{gs_items[row].DistinctPeptideCount}</Cell>}
                />
                <Column
                  name="Peptides"
                  cellRenderer={row => <Cell key={row+''}>{gs_items[row].Peptides}</Cell>}
                />
              </Table>
            )
          })
    .build()

export const GlycoSightFileURLNode = MetaNode("GlycoSightURL")
    .meta({
        label:"GlycoSight File URL",
        description: "URL to storage of MZID archives",
        icon: [glygen_icon],
        external: true,
    })
    .codec(GlycoSightFileURLC)
    .view((data) => (
        <div>
            <h2 className="prose max-w-none">GlycoSight File: {data.url}</h2>
        </div>
    ))
    .build()

export const GlycoSightProcessNode = MetaNode("GlycoSightProcessNode")
    .meta({
        label: `Launch GlycoSight Analysis`,
        description: `Run GlycoSight Analysis on MZID file(s)`,
        icon: [glygen_icon],
        external: true,
    })
    .inputs({ file: GlycoSightFileURLNode })
    .output(GlycoSightOutputNode)
    .resolve( async ( props ) => {
        
      const fileStream = await fileAsStream(props.inputs.file);
      const chunks = [];
      for await (const chunk of fileStream) {
        chunks.push(Buffer.from(chunk));
      }
      const fileBuffer = Buffer.concat(chunks);

      const result = await UploadAndAnalyze(props.inputs.file, fileBuffer);
      
      return result;
    })
    .story((props) => ({
      abstract: `The N-Linked glycan sites analysis was launched through the GlycoSight servers.`,
    }))
    .build()

export const GlycoSightFileUpload = MetaNode("GlycoSightUpload")
    .meta({
        label: "Upload GlycoSight MZID",
        description: "Upload MZID data for peptide-to-Uniprot sample mapping",
        icon: [glygen_icon],
        external: true,
    })
    .codec(FileC)
    .inputs() 
    .output(GlycoSightFileURLNode)
    .prompt(props => <>
        {/* XXX: <FilePrompt {...props} /> */}
        TODO: <FilePrompt {...props} example={clientLoadExample} />
        {props.output ? 
            <SafeRender 
                component={FileURL.view} 
                props={props.output} 
            /> : 
            null}</>
        )
        .resolve(async (props) => props.data)
        .story((props) => ({
          abstract: `MZID glycosylation data for N-linked glycan analysis was uploaded.`,
          introduction: `GlycoSight parses peptide identification results, in mzIdentML format, to find putative N-glycosylation sites. The analyzed sample must have been prepared with PNGase-F, an enzyme which cleaves N-glycans from the Asn amino-acid and deamidates it, and changes the molecular weight of the Asn residue by +0.98. The output table provides the UniProt accession, the amino-acid position of the Asn exibiting deamidation, the number of tandem mass spectra which show deamidation of the residue, the number of distinct peptide sequences which show deamindation of the residue, and the peptide sequences themselves. Peptide sequences are separated by a comma, and the corresponding Asn residue shown in lowercase.`,
          methods: `First, peptides were identified and mapped to UniProt Accession ID numbers. The peptides were then analyzed for the mass shift of +0.98, indicative of deglycosylation following treatment with the enzyme PNGase-F.`,
          legend: `A table listing identified peptide signatures, or a comma-separated list of peptide signatures, UniProt Protein Accession IDs, associated amino acid position, spectral count, and distinct peptide counts for each analyzed peptide signature. Putative N-glycosylation sites are annotated with lower-case "n" to represent the deamidated asparagine (Asn) residue.`,
        }))
        .build()




// TODO: To pick up at a later date
//     .resolve(async (props) => {
//         const result = resolvePython(props.data, props.notify);
//         return result;
//         }
//     )
//     .story(props => "We did it!")
//     .build()

// async function resolvePython(data, notify) {
//     await python(
//         'components.gly_gen.glycosight.load_mzid_file', 
//         { kargs: [data] },
//         message => notify({ type: 'info', message }),
//     )
// }
