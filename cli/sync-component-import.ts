import fs from 'fs'
import path from 'path'
fs.writeFileSync('components/index.ts',
  fs.readdirSync(path.join(__dirname, '..', 'components'))
    .map(p => path.basename(p))
    .filter(component => component !== 'index.ts')
    .map(component => `export * from ${JSON.stringify('./' + path.basename(component))}`)
    .join('\n')
)
