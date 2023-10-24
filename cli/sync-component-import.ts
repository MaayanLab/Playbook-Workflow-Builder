import fs from 'fs'
import glob from 'glob'
import path from 'path'

const base = path.join(__dirname, '..', 'components')
const components = glob.sync(path.join(base, '**', 'package.json').split(path.sep).join(path.posix.sep))
  .filter(p => !p.includes('node_modules'))
  .map(p => path.dirname(p))
components.sort()

fs.writeFileSync(path.join(base, 'index.ts'), [
  `import { MetaNode, MetaNodesFromExports } from "@/spec/metanode"`,
  `export const components: string[] = []`,
  `export const metanodes: MetaNode[] = []`,
  ...components
    .flatMap(component => {
      const componentPath = path.relative(base, component)
      const componentSlug = componentPath.replaceAll(/[/-]/g, '_')
      return [
        `import * as ${componentSlug} from ${JSON.stringify(`./${componentPath}`)}`,
        `import ${componentSlug}_packageJson from ${JSON.stringify(`./${componentPath}/package.json`)}`,
        `metanodes.push(...MetaNodesFromExports(${componentSlug}, ${componentSlug}_packageJson))`,
        `components.push(${JSON.stringify(componentPath)})`,
      ]
    }),
  ].join('\n')
)
