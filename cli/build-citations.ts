import fs from 'fs'
import path from 'path'
import process from 'process'
import * as glob from 'glob'
import * as dict from '@/utils/dict'
import * as array from '@/utils/array'
import currentCitations from '@/utils/citations/citations'

function re_findall(RE: RegExp, doc: string) {
  const re = new RegExp(RE)
  const M = [] as RegExpExecArray[]
  let m
  while ((m = re.exec(doc)) !== null) {
    M.push(m)
  }
  return M
}

function locateCitations() {
  return ['*.ts', '*.tsx'].flatMap(ext =>
    [...glob.sync(path.join(__dirname, '..', 'components', '**', ext).split(path.sep).join(path.posix.sep))]
  ).flatMap(p =>
    re_findall(/doi:([0-9a-z\./_-]+)/gi, fs.readFileSync(p).toString())
  ).map(m => m[1])
}

async function fetchCitation(doi: string) {
  console.info(`Fetching ${doi}...`)
  const req = await fetch(`https://citation.crosscite.org/format?style=nature&lang=en-US&doi=${encodeURIComponent(doi)}`, {
    method: 'GET',
    headers: {
      Accept: 'text/plain',
    },
  })
  const res = await req.text()
  return res.replace(/^1\.\s*/, '').replace(/\s*$/g, '')
}

async function writeCitations(citations: { [id: string]: string }) {
  if (dict.isEmpty(citations)) {
    console.info(`Nothing to do.`)
    return
  }
  console.info(`Adding ${dict.keys(citations).join(', ')}...`)
  Object.assign(currentCitations, citations)
  await new Promise<void>((resolve, reject) =>
    fs.writeFile(path.join(__dirname, '..', 'utils', 'citations', 'citations.ts'), [
      'export default {',
      ...dict.sortedItems(currentCitations).map(({ key, value }) => `  ${JSON.stringify(key)}: ${JSON.stringify(value)},`),
      '}',
    ].join('\n'), (err) => {
      if (err) reject(err)
      else resolve()
    })
  )
}

Promise.all(
  array.unique(
    locateCitations().filter(doi => !(doi in currentCitations))
  ).map(doi =>
    fetchCitation(doi).then(citation => ({ key: doi, value: citation }))
  )
).then(citations =>
  dict.init(citations)
).then(citations =>
  writeCitations(citations)
).catch((error) => {
  console.error(error)
  process.exit(1)
})
