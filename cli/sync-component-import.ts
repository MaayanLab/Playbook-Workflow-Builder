import fs from 'fs'
import glob from 'glob'
import path from 'path'

const base = path.join(__dirname, '..', 'components')
const components = glob.sync(path.join(base, '**/package.json'))
  .map(p => path.dirname(p))
components.sort()

fs.writeFileSync(path.join(base, 'index.ts'),
  components
    .map(component => `export * from ${JSON.stringify(`./${path.relative(base, component)}`)}`)
    .join('\n')
)
