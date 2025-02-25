/**
 * This converts PWB workflows and components into Research Object Crates
 */
import type { FPL } from "@/core/FPPRG"
import type KRG from "@/core/KRG"
import * as dict from '@/utils/dict'
import packageJson from '@/package.json'
import jsonld from 'jsonld'
import YAML from 'yaml'
import { Author, Metadata, fpl_expand, toISO8601TimeString } from "./common";
import { cwl_for_playbook } from "./fpl2cwl"
import FPL2BCO from "./fpl2bco"
import fpl2json from "./fpl2json"

const version = packageJson.version

export async function fpl2ro_crate_metadata(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null, cwl_files?: Record<string, any> }) {
  const { processLookup, story } = await fpl_expand(props)
  let doc: string
  if (props.metadata?.description) {
    doc = props.metadata.description
  } else {
    doc = story.ast.flatMap(part => part.tags.includes('abstract') ? [part.type === 'bibitem' ? '\n' : '', part.text] : []).join('')
  }
  return jsonld.flatten({
    "@context": "https://w3id.org/ro/crate/1.1/context",
    "@graph": [
      {
        "@id": "ro-crate-metadata.json",
        "@type": "CreativeWork",
        "about": {
          "@id": "./"
        },
        "conformsTo": {
          "@id": "https://w3id.org/ro/crate/1.1"
        }
      },
      {
        "@id": "./",
        "@type": "Dataset",
        "conformsTo": [
          {
            "@id": "https://w3id.org/ro/crate/1.1"
          },
          {
            "@id": "https://w3id.org/workflowhub/workflow-ro-crate/1.0"
          }
        ],
        "license": 'https://creativecommons.org/licenses/by/4.0/',
        "name": props.metadata?.title ?? 'Playbook',
        "description": doc,
        "isBasedOn": `${process.env.PUBLIC_URL}/report/${props.fpl.id}`,
        "mainEntity": {
          "@id": "playbook.json",
        },
        "hasPart": [
          {"@id": "playbook.json"},
          {"@id": "workflow.cwl"},
          {"@id": "bco.json"},
        ],
      },
      {
        "@id": 'playbook.json',
        "@type": [
          "File",
          "SoftwareSourceCode",
          "ComputationalWorkflow",
        ],
        "dct:conformsTo": "https://bioschemas.org/profiles/ComputationalWorkflow/1.0-RELEASE/",
        "url": `${process.env.PUBLIC_URL}/report/${props.fpl.id}`,
        "encodingFormat":"application/json",
        "contentUrl": `${process.env.PUBLIC_URL}/api/db/fpl/${props.fpl.id}/output/export`,
        "programmingLanguage": {
          "@id": "#playbook_workflow_builder_workflow",
          "@type": "ComputerLanguage",
          "name": "Playbook Workflow Builder Workflow",
          "url": {
            "@id": "https://playbook-workflow-builder.cloud"
          },
          "version": version,
          "alternateName": "PWB",
        },
        "license": 'https://creativecommons.org/licenses/by/4.0/',
        "dateCreated": toISO8601TimeString(),
        "version": 0,
        "sdPublisher": {
          "@id": "#playbook_workflow_builder_cloud",
          "@type": "Organization",
          "url": "https://playbook-workflow-builder.cloud"
        },
        "creator": props.author?.orcid ? [
          {
            "@id": `#${props.author.orcid}`,
            "@type": "Person",
            "name": props.author.name,
            // "email": props.author.email,
            "affiliation": props.author.affiliation ? {
              "@type": "Organization",
              "name": props.author.affiliation
            } : undefined,
          },
        ] : undefined,
        "input": dict.values(processLookup).filter(({ metanode }) => 'codec' in metanode).map(({ index, node, metanode }) => ({
          "@type": "FormalParameter",
          "dct:conformsTo": "https://bioschemas.org/profiles/FormalParameter/1.0-RELEASE/",
          "name": metanode.meta.label,
          "description": metanode.meta.description,
          "encodingFormat": "application/json",
        })),
        "output": dict.values(processLookup).map(({ index, node, metanode }) => ({
          "@type": "FormalParameter",
          "dct:conformsTo": "https://bioschemas.org/profiles/FormalParameter/1.0-RELEASE/",
          "name": metanode.output.meta.label,
          "description": metanode.output.meta.description,
          "encodingFormat": "application/json",
        })),
      },
      {
        "@id": "bco.json",
        "@type": [
          "File",
          "SoftwareSourceCode"
        ],
        "programmingLanguage": {
          "@type": "ComputerLanguage",
          "name": "BioCompute Object",
          "alternateName": "BCO",
          "identifier": "https://w3id.org/ieee/ieee-2791-schema",
          "url": "https://biocomputeobject.org/",
        },
        "encodingFormat":"application/json",
        "contentUrl": `${process.env.PUBLIC_URL}/api/bco/${props.fpl.id}`,
        "name": "Playbook BCO",
      },
      {
        "@id": "#cwl",
        "@type": "ComputerLanguage",
        "name": "Common Workflow Language",
        "alternateName": "CWL",
        "identifier": "https://w3id.org/cwl/v1.0/",
        "url": "https://www.commonwl.org/",
      },
      ...(props.cwl_files ? dict.items(props.cwl_files ?? {}).map(({ key, value }) =>
        key.endsWith('.cwl') ? {
          "@id": key,
          "@type": [
            "File",
            "SoftwareSourceCode"
          ],
          "programmingLanguage": {
            "@id": "#cwl"
          },
          "contentUrl": `${process.env.PUBLIC_URL}/api/v1/cwl/${props.fpl.id}/${key}`,
          "contentSize": value.length,
        }
        : key.endsWith('.yaml') ? {
          "@id": key,
          "@type": [
            "File",
            "SoftwareSourceCode"
          ],
          "programmingLanguage": {
            "@id": "#yaml"
          },
          "contentUrl": `${process.env.PUBLIC_URL}/api/v1/cwl/${props.fpl.id}/${key}`,
          "contentSize": value.length,
        }
        : {
          "@id": key,
          "@type": "File",
          "contentUrl": `${process.env.PUBLIC_URL}/api/v1/cwl/${props.fpl.id}/${key}`,
          "contentSize": value.length,
        }
      ) : [
        {
          "@id": 'workflow.cwl',
          "@type": [
            "File",
            "SoftwareSourceCode"
          ],
          "programmingLanguage": {
            "@id": "#cwl"
          },
          "contentUrl": `${process.env.PUBLIC_URL}/api/v1/cwl/${props.fpl.id}/workflow.cwl`,
        },
        {
          "@id": 'inputs.yaml',
          "@type": [
            "File",
            "SoftwareSourceCode"
          ],
          "programmingLanguage": {
            "@id": "#yaml"
          },
          "contentUrl": `${process.env.PUBLIC_URL}/api/v1/cwl/${props.fpl.id}/inputs.yaml`,
        },
      ]),
    ],
  } as any)
}

export async function fpl2ro_crate(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null }) {
  const cwl_content = await cwl_for_playbook(props)
  const cwl_files = dict.init(
    dict.items(cwl_content)
      .map(({ key, value }) => ({ key, value: YAML.stringify(value) }))
  )
  return {
    ...cwl_files,
    'playbook.json': JSON.stringify(await fpl2json(props)),
    'bco.json': JSON.stringify(await FPL2BCO(props)),
    'ro-crate-metadata.json': JSON.stringify(await fpl2ro_crate_metadata({ ...props, cwl_files })),
  }
}
