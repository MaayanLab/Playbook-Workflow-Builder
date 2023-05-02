import { Configuration, OpenAIApi } from "openai"
import KRG from '@/core/KRG'
import * as components from '@/components'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { zodToTs, createTypeAlias, printNode } from 'zod-to-ts'
import { MetaNode, ProcessMetaNode } from '@/spec/metanode'

/**
 * Load MetaNodes
 */
const krg = new KRG()
dict.values(components as any).flatMap(component => array.ensureArray(component))
  .filter((component): component is MetaNode => typeof component === 'object' && 'spec' in component)
  .forEach(metanode => krg.add(metanode))

/**
 * Initialize OpenAI API
 */
const configuration = new Configuration({
  apiKey: process.env.OPENAI_API_KEY,
})
const openai = new OpenAIApi(configuration)

function slugify(spec: string) {
  return spec.replaceAll(/[\[\]\(\)\{\} ,-]/g, '_').replaceAll(/_+$/g, '')
}

function try_story<T extends ProcessMetaNode>(metanode: T, props: T['story'] extends (undefined | ((arg: infer A) => string)) ? A : never) {
  try {
    if (metanode.story) {
      return metanode.story(props as any)
    }
  } catch (e) {}
  return ''
}

/**
 * Create typescript types for each Data MetaNode
 */
const DataNodeSpecs = dict.init(
  krg.getDataNodes()
    .map(metanode => {
      const typename = slugify(metanode.spec)
      const { node } = zodToTs(metanode.zod, typename)
      const typealias = createTypeAlias(node, typename)
      return {
        key: metanode.spec,
        value: [
          `/**`,
          ` * ${metanode.meta.label}: ${metanode.meta.description}`,
          ` */`,
          printNode(typealias),
          '',
        ].join('\n')
      }
    })
)

/**
 * Create typescript functions for each Process MetaNode
 */
const ProcessNodeSpecs = dict.init(
  krg.getProcessNodes()
    .map(metanode => {
      return {
        key: metanode.spec,
        value: [
          ...array.unique(
            dict.values(metanode.inputs)
              .flatMap(value => array.ensureArray(value))
              .map(value => value.spec)
          ).map(spec => DataNodeSpecs[spec]),
          `/**`,
          ` * ${metanode.meta.label}: ${metanode.meta.description}`,
          ` * ${try_story(metanode, {
                inputs: dict.init(
                  dict.items(metanode.inputs).map(({ key, value }) => ({
                    key, value: `[${Array.isArray(value) ? `${slugify(value[0].spec)}[]` : slugify(value.spec)}]`,
                  }))
                ),
                output: `[${slugify(metanode.output.spec)}]`,
              })}`,
          ` */`,
          `async function ${slugify(metanode.spec)}(${
            dict.items(metanode.inputs)
              .map(({ key, value }) =>
                `${slugify(key)}: ${Array.isArray(value) ? `${slugify(value[0].spec)}[]` : slugify(value.spec)}`
              ).join(', ')
          }): Promise<${slugify(metanode.output.spec)}>`,
          '',
        ].join('\n')
      }
    })
)

/**
 * Create more elaborate context for a metanode, including the node itself and its next processes
 */
const S2ProcessNodeSpecFactory = () => dict.init(
  krg.getProcessNodes()
    .map(metanode => {
      const nextProcessNodes = krg.getNextProcess(metanode.output.spec)
      return {
        key: metanode.spec,
        value: [
          ...array.unique(
            [metanode, ...nextProcessNodes].flatMap(metanode =>
              dict.values(metanode.inputs)
                .flatMap(value => array.ensureArray(value))
                .map(value => value.spec)
            )
          ).map(spec => DataNodeSpecs[spec]),
          ...[metanode, ...nextProcessNodes].map(processNode => ProcessNodeSpecs[processNode.spec]),
          `async function main({`,
          ...dict.keys(metanode.inputs).map(key =>
            `  ${slugify(key)},`
          ),
          `}: {`,
          ...dict.items(metanode.inputs).map(({ key, value }) =>
            `  ${slugify(key)}: ${Array.isArray(value) ? `${slugify(value[0].spec)}[]` : slugify(value.spec)},`,
          ),
          `}) {`,
          `  const x: ${slugify(metanode.output.spec)} = await ${slugify(metanode.spec)}({${dict.keys(metanode.inputs).map(key => slugify(key)).join(', ')}})`,
          ...nextProcessNodes.map(nextMetanode =>
            `  console.log(await ${slugify(nextMetanode.spec)}({${
              dict.items(nextMetanode.inputs).map(({ key, value }) =>
                metanode.output.spec === array.ensureArray(value)[0].spec ?
                  `${slugify(key)}: ${Array.isArray(value) ? `[x]` : `x`}`
                  : slugify(key)
              ).join(', ')}}))`
          ),
          `}`,
        ].join('\n')
      }
    })
)

async function main(dry_run = true, s2 = false) {
  const Specs = s2 ? S2ProcessNodeSpecFactory() : ProcessNodeSpecs
  const keys = dict.keys(Specs)
  console.log(['id','embedding'].join('\t'))
  if (dry_run) {
    keys.forEach(key => {
      console.log([key, JSON.stringify(Specs[key])].join('\t'))
    })
  } else {
    const embeddings = await openai.createEmbedding({
      model: 'text-embedding-ada-002',
      input: keys.map(key => Specs[key]),
    })
    embeddings.data.data.forEach((datum, i) => {
      console.log([keys[datum.index], JSON.stringify(datum.embedding)].join('\t'))
    })
  }
}
main(
  !process.argv.slice(2).some(arg => arg === '--embed'),
  process.argv.slice(2).some(arg => arg === '--s2'),
)
