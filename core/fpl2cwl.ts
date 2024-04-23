/**
 * This converts PWB workflows and components into CWL workflows, arguments & CommandLineTool specifications.
 */
import type { ProcessMetaNode } from "@/spec/metanode";
import type { FPL } from "@/core/FPPRG"
import type KRG from "@/core/KRG"
import * as array from '@/utils/array'
import * as dict from '@/utils/dict'
import packageJson from '@/package.json'
import { Author, Metadata, fpl_expand } from "./common";

const docker_tag = 'maayanlab/playbook-partnership'
const version = packageJson.version

export function cwl_cli_for_component(component: ProcessMetaNode) {
  return {
    cwlVersion: 'v1.0',
    class: 'CommandLineTool',
    label: component.meta.label,
    doc: component.meta.description,
    baseCommand: [
      '/app/cli/pwb.sh',
      component.spec,
    ],
    requirements: [
      {
        class: 'DockerRequirement',
        dockerImageId: `${docker_tag}:${version}`,
      },
    ],
    inputs: {
      ...('codec' in component ? {
        data: {
          type: 'File',
          inputBinding: {
            prefix: `--data=`,
            separate: false,
          },
        },
      } : {}),
      ...dict.init(dict.items(component.inputs).map(({ key, value }) => ({
        key: `inputs.${key}`,
        value: {
          label: array.ensureArray(value)[0].meta.label,
          doc: array.ensureArray(value)[0].meta.description,
          type: Array.isArray(value) ? 'File[]' : 'File',
          inputBinding: {
            prefix: `--inputs.${key}=`,
            separate: false,
          },
        },
      }))),
      outputFilename: {
        type: 'string',
        default: 'output.json',
        inputBinding: {
          prefix: `--output=`,
          separate: false,
        },
      },
    },
    outputs: {
      output: {
        label: component.output.meta.label,
        doc: component.output.meta.description,
        type: 'File',
        outputBinding: {
          glob: '$(inputs.outputFilename)',
        },
      },
    },
  }
}

export async function cwl_for_playbook(props: { krg: KRG, fpl: FPL, metadata?: Metadata, author?: Author | null }) {
  const { processLookup, story } = await fpl_expand(props)
  return {
    ...dict.init(array.unique(dict.values(processLookup).map(({ metanode }) => metanode)).map(metanode => ({
      key: `${metanode.spec}.cwl`,
      value: cwl_cli_for_component(metanode),
    }))),
    'workflow.cwl': {
      cwlVersion: 'v1.2',
      class: 'Workflow',
      label: props.metadata?.title,
      doc: props.metadata?.description ?? story,
      requirements: {},
      inputs: dict.init(dict.values(processLookup).filter(({ metanode }) => 'codec' in metanode).map(({ index, node, metanode }) => ({
        key: `step-${index+1}-data`,
        value: {
          label: metanode.meta.label,
          doc: metanode.meta.description,
          type: 'File',
        },
      }))),
      outputs: dict.init(dict.values(processLookup).map(({ index, node, metanode }) => ({
        key: `step-${index+1}-output`,
        value: {
          label: metanode.output.meta.label,
          doc: metanode.output.meta.description,
          type: 'File',
          outputSource: `step-${index+1}/output`,
        },
      }))),
      steps: dict.init(dict.values(processLookup).map(({ index, node, metanode }) => ({
        key: `step-${index+1}`, value: {
          run: `${metanode.spec}.cwl`,
          label: metanode.meta.label,
          doc: metanode.meta.description,
          in: {
            ...('codec' in metanode ? {
              data: {
                source: `step-${index+1}-data`
              },
            } : {}),
            ...dict.init(dict.items(node.inputs).map(({ key, value }) => ({ key, value: processLookup[value.id] })).map(({ key, value }) => ({
              key: `inputs.${key}`,
              value: {
                label: value.metanode.output.meta.label,
                source: `step-${value.index+1}/output`,
              },
            }))),
            outputFilename: {
              default: `step-${index+1}-output.json`,
            },
          },
          out: ['output'],
        },
      }))),
    },
    'inputs.yaml': dict.init(await Promise.all(dict.values(processLookup).filter(({ metanode }) => 'codec' in metanode).map(async ({ index, node, metanode }) => ({
      key: `step-${index+1}-data`,
      value: {
        class: 'File',
        contents: JSON.stringify((await node.output())?.value ?? null),
      },
    }))))
  }
}
