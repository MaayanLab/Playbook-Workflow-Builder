import KRG from '@/core/KRG'
import { metanodes } from '@/components'
import * as dict from '@/utils/dict'

const krg = new KRG()
metanodes.forEach(metanode => krg.add(metanode))

console.log(
  [
    ['id', 'type', 'input', 'output', 'label', 'description', 'color', 'abstract', 'legend', 'introduction', 'methods', 'results'].join('\t'),
    ...[
      ...krg.getDataNodes(),
      ...krg.getPromptNodes(),
      ...krg.getResolveNodes(),
    ].map(node => {
      const story = 'story' in node ? node.story({}) : {}
      return [
        node.spec,
        'prompt' in node ? 'Prompt'
          : node.kind === 'process' ? 'Resolver'
          : 'Data',
        node.kind === 'process' ? JSON.stringify(dict.init(dict.items(node.inputs).map(({ key, value }) => ({ key, value: Array.isArray(value) ? [value[0].spec] : value.spec })))) : '',
        node.kind === 'process' ? node.output.spec : '',
        node.meta.label,
        node.meta.description,
        node.meta.color,
        story.abstract ?? '',
        story.legend ?? story.figureLegend ?? story.tableLegend ?? '',
        story.introduction ?? '',
        story.methods ?? '',
        story.results ?? '',
      ].join('\t')
    }),
  ].join('\n')
)
