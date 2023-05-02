import fs from 'fs'
import glob from 'glob'
import path from 'path'

const base = path.join(__dirname, '..', 'components')
const components = glob.sync(path.join(base, '**', 'package.json').replace(/\\/g, '/'))
  .filter(p => !p.includes('node_modules'))
  .map(p => path.dirname(p))
components.sort()

fs.writeFileSync(path.join(base, 'index.ts'), [
  `export const components: string[] = []`,
  ...components
    .flatMap(component => {
      const componentPath = path.relative(base, component)
      return [
        `export * from ${JSON.stringify(`./${componentPath}`)}`,
        `components.push(${JSON.stringify(componentPath)})`,
      ]
    }),
  ].join('\n')
)
