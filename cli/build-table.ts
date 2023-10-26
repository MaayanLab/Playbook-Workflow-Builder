import KRG from '@/core/KRG'
import * as components from '@/components'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import { MetaNode } from '@/spec/metanode'

const krg = new KRG()
dict.values(components as any).flatMap(component => array.ensureArray(component))
  .filter((component): component is MetaNode => typeof component === 'object' && 'spec' in component)
  .forEach(metanode => krg.add(metanode))

console.log(
  [
    ['id', 'type', 'input', 'output', 'label', 'description', 'color'].join('\t'),
    ...[
      ...krg.getDataNodes(),
      ...krg.getPromptNodes(),
      ...krg.getResolveNodes(),
    ].map(node => [
      node.spec,
      'prompt' in node ? 'Prompt'
        : node.kind === 'process' ? 'Resolver'
        : 'Data',
      node.kind === 'process' ? JSON.stringify(dict.init(dict.items(node.inputs).map(({ key, value }) => ({ key, value: Array.isArray(value) ? [value[0].spec] : value.spec })))) : '',
      node.kind === 'process' ? node.output.spec : '',
      node.meta.label,
      node.meta.description,
      node.meta.color,
    ].join('\t')),
  ].join('\n')
)
