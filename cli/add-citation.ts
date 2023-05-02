import fs from 'fs'
import path from 'path'
import process from 'process'
import * as dict from '@/utils/dict'
import citations from '@/utils/citations/citations'

async function fetchCitation(doi: string) {
  const req = await fetch(`https://citation.crosscite.org/format?style=nature&lang=en-US&doi=${encodeURIComponent(doi)}`, {
    method: 'GET',
    headers: {
      Accept: 'text/plain',
    },
  })
  const res = await req.text()
  return res.replace(/^1\./, '').replace(/\n*$/g, '')
}

const [_process, _script, doi] = process.argv

if (!(doi in citations)) {
  fetchCitation(doi)
    .then((citation) => {
      console.info(`Adding ${doi}...`)
      Object.assign(citations, { [doi]: citation })
      fs.writeFileSync(path.join(__dirname, '..', 'utils', 'citations', 'citations.ts'), [
        'export default {',
        ...dict.sortedItems(citations).map(({ key, value }) => `  ${JSON.stringify(key)}: ${JSON.stringify(value)},`),
        '}',
      ].join('\n'))
    })
    .catch((error) => {
      console.error(error)
      process.exit(1)
    })
} else {
  console.warn(`${doi} already registered`)
}
