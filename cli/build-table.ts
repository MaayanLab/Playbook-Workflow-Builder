import KRG from '@/core/KRG'
import { metanodes } from '@/components'
import * as dict from '@/utils/dict'
import { ProcessMetaNode } from '@/spec/metanode'

const krg = new KRG()
metanodes.forEach(metanode => krg.add(metanode))

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

console.log(
  [
    ['id', 'slug', 'type', 'func', 'label', 'description', 'story'].join('\t'),
    ...[
      ...krg.getDataNodes(),
      ...krg.getPromptNodes(),
      ...krg.getResolveNodes(),
    ].map(node => [
      node.spec,
      slugify(node.spec),
      'prompt' in node ? 'Prompt'
        : node.kind === 'process' ? 'Resolver'
        : 'Data',
        node.kind === 'process' ? `${slugify(node.spec)}(${[
          ...dict.items(node.inputs).map(({ key, value }) => `${key}: ${Array.isArray(value) ? `list[${slugify(value[0].spec)}]` : slugify(value.spec)}`),
          ...('prompt' in node ? ['input: str'] : []),
        ].join(', ')}) -> ${slugify(node.output.spec)}` : '',
      node.meta.label,
      node.meta.description,
      'story' in node ? try_story(node, {
        inputs: dict.init(
          dict.items(node.inputs).map(({ key, value }) => ({
            key, value: `{${Array.isArray(value) ? `[${slugify(value[0].spec)}]` : slugify(value.spec)}}`,
          }))
        ),
        output: `{${slugify(node.output.spec)}}`,
      }) : '',
    ].join('\t')),
  ].join('\n')
)
