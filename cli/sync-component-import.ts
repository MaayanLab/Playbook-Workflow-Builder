import fs from 'fs'
import path from 'path'

const components = fs.readdirSync(path.join(__dirname, '..', 'components'))
  .filter(p => !p.startsWith('_') && fs.lstatSync(path.join(__dirname, '..', 'components', p)).isDirectory())
  .map(p => path.basename(p))
components.sort()

fs.writeFileSync('components/index.ts',
  components
    .map(component => `export * from ${JSON.stringify('./' + path.basename(component))}`)
    .join('\n')
)
