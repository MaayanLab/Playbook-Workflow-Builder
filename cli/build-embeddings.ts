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

async function main() {
  const keys = dict.keys(ProcessNodeSpecs)
  console.log(['id','embedding'].join('\t'))
  const embeddings = await openai.createEmbedding({
    model: 'text-embedding-ada-002',
    input: keys.map(key => ProcessNodeSpecs[key]),
  })
  embeddings.data.data.forEach((datum, i) => {
    console.log([keys[datum.index], JSON.stringify(datum.embedding)].join('\t'))
  })
}
main()
