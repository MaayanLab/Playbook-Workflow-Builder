import { z } from "zod"
import type KRG from "@/core/KRG"
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { FPL } from "./FPPRG"

const uriSchema = z.object({
  uri: z.string(),
  access_time: z.string().optional(),
  filename: z.string().optional(),
})
export const BCOschema = z.object({
  etag: z.string(),
  object_id: z.string(),
  spec_version: z.string(),
  usability_domain: z.array(z.string()),
  provenance_domain: z.object({
    embargo: z.object({}),
    name: z.string(),
    version: z.string(),
    license: z.string(),
    derived_from: z.string(),
    contributors: z.array(
      z.object({
        name: z.string(),
        orcid: z.string().optional(),
        affiliation: z.string(),
        contribution: z.array(z.string()),
        email: z.string()
      }),
    ),
    review: z.array(
      z.object({
        reviewer: z.object({
          name: z.string(),
          orcid: z.string(),
          affiliation: z.string(),
          contribution: z.array(z.string()),
          email: z.string()
        }),
        status: z.string()
      })
    ),
    created: z.string(),
    modified: z.string()
  }),
  description_domain: z.object({
    keywords: z.array(z.string()),
    platform: z.array(z.string()),
    pipeline_steps: z.array(
      z.object({
        name: z.string(),
        description: z.string().optional(),
        version: z.string().optional(),
        step_number: z.number(),
        input_list: z.array(uriSchema),
        output_list: z.array(uriSchema),
        prerequisite: z.array(
          z.object({
            name: z.string(),
            uri: uriSchema,
          })
        ).optional()
      }),
    )
  }),
  parametric_domain: z.array(
    z.object({ step: z.string(), param: z.string(), value: z.string() })
  ),
  execution_domain: z.object({
    external_data_endpoints: z.array(
      z.object({ name: z.string(), url: z.string() })
    ),
    software_prerequisites: z.array(
      z.object({
        name: z.string(),
        version: z.string(),
        uri: uriSchema,
      }),
    ),
    environment_variables: z.record(z.string(), z.string()),
    script_driver: z.string(),
    script: z.array(
      z.object({ uri: uriSchema })
    )
  }),
  io_domain: z.object({
    input_subdomain: z.array(z.object({ uri: uriSchema, mediatype: z.string().optional() })),
    output_subdomain: z.array(z.object({ uri: uriSchema, mediatype: z.string().optional() })),
  }),
  error_domain: z.object({
    algorithmic_error: z.object({}).optional(),
    empirical_error: z.object({}).optional(),
  })
})

function toBCOTimeString(date?: Date) {
  if (date === undefined) date = new Date()
  return date.toISOString().replace(/Z$/, '000')
}
type PromiseType<T> = T extends Promise<infer RT> ? RT : never
type FullFPL = Array<PromiseType<ReturnType<FPL['toJSONWithOutput']>>>

export default function FPL2BCO(krg: KRG, fpl: FullFPL): z.infer<typeof BCOschema> {
  const processLookup = dict.init(fpl.map((step, index) => ({
    key: step.process.id,
    value: {
      index,
      node: step.process,
      metanode: krg.getProcessNode(step.process.type)
    }
  })))
  return {
    etag: '',// TODO
    object_id: '', // TODO
    spec_version: 'https://w3id.org/ieee/ieee-2791-schema/',
    usability_domain: [
      // TODO
      'Some description about this workflow',
    ],
    provenance_domain: {
      embargo: {}, // ?
      name: 'Playbook Partnership',
      version: '1.0',
      license: 'https://creativecommons.org/licenses/by-nc-sa/4.0/',
      derived_from: 'NA',
      contributors: [
        {
          name: 'Daniel Clarke',
          orcid: 'https://orcid.org/0000-0003-3471-7416',
          affiliation: 'Mount Sinai',
          contribution: [
            'authoredBy',
          ],
          email: 'danieljbclarkemssm@gmail.com',
        },
      ],
      review: [],
      created: toBCOTimeString(), // TODO: datetime
      modified: toBCOTimeString(), // TODO: datetime
    },
    description_domain: {
      keywords: array.unique(
        Object.values(processLookup)
          .flatMap(({ metanode }) =>
            metanode.meta.tags ? Object.keys(metanode.meta.tags) : []
          )
      ),
      platform: ['Debian GNU/Linux 11'],
      pipeline_steps: Object.values(processLookup).map(({ index, node, metanode }) => ({
        name: metanode.meta.label,
        description: metanode.meta.description,
        // version: metanode.meta.version,
        step_number: index + 1,
        prerequisite: Object.values(node.inputs).map(input => processLookup[input.id]).map(inputProcess => ({
          name: `Output of step ${inputProcess.index+1}`,
          uri: {
            uri: `#/${inputProcess.index}/process/output`,
            access_time: toBCOTimeString(),
          },
        })),
        input_list: Object.values(node.inputs).map(input => processLookup[input.id]).map(inputProcess => ({
          uri: `#/${inputProcess.index}/process/output`,
          access_time: toBCOTimeString(),
        })),
        output_list: [
          {
            uri: `#/${index}/process/output`,
            access_time: toBCOTimeString(),
          }
        ],
      })),
    },
    parametric_domain: Object.values(processLookup).filter(({ node }) => node.data !== null).map(({ index, node, metanode }) => ({
      step: `${index+1}`,
      param: 'stdin',
      value: node.data ? node.data.value : '',
    })),
    execution_domain: {
      // TODO: load prereqs from steps in use (?)
      external_data_endpoints: [],
      software_prerequisites: [
        {
          name: 'Docker',
          version: '20.10.21',
          uri: {
            access_time: toBCOTimeString(),
            uri: 'https://docs.docker.com/get-docker/',
          },
        },
      ],
      environment_variables: {},
      script_driver: 'shell',
      script: [
        {
          uri: {
            uri: 'https://github.com/nih-cfde/playbook-partnership/tree/dev/cli/playbook-partnership-executor.sh',
            filename: 'playbook-partnership-executor.sh',
          }
        }
      ]
    },
    io_domain: {
      input_subdomain: [
        {
          uri: {
            uri: `${process.env.PUBLIC_URL}/api/db/fpl/${fpl[fpl.length-1].id}`,
          },
          mediatype: 'application/json',
        },
      ],
      output_subdomain: [
        {
          uri: {
            uri: `${process.env.PUBLIC_URL}/api/db/fpl/${fpl[fpl.length-1].id}/output`,
          },
          mediatype: 'application/json'
        }
      ]
    },
    error_domain: {
      // TODO
      algorithmic_error: {},
      empirical_error: {},
    },
  }
}
