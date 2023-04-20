import fs from 'fs'

export async function download(file: { url: string }) {
  if (file.url.startsWith('file://')) {
    return fs.readFileSync(file.url.slice('file://'.length)).toString()
  } else {
    const req = await fetch(file.url)
    return await req.text()
  }
}
