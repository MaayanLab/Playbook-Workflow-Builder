import { spawn } from 'child_process'
import fs from 'fs/promises'
import os from 'os'
import path from 'path'
import JSZip from 'jszip'
import fpprg from '@/app/fpprg'
import { getKRG } from '@/app/server-krg'
import FPL2TEX from '@/core/fpl2tex'
import { getRequestSessionWithId } from '@/app/extensions/auth/server'
import { getUser, getUserOrcidAccount } from '@/app/extensions/auth/user'
import { NotFoundError } from '@/spec/error'
import { routeResponse } from '@/app/app/api/_lib/response'

function latexmk(main: string, cwd: string) {
  return new Promise<void>((resolve, reject) => {
    const proc = spawn('latexmk', ['-pdf', '-cd', '-lualatex', '-interaction=nonstopmode', '-f', main], { cwd })
    let stderr = ''
    proc.stdout.on('data', data => process.stdout.write(`[latexmk] ${data}`))
    proc.stderr.on('data', data => { stderr += data; process.stderr.write(`[latexmk] ${data}`) })
    proc.on('close', code => code === 0 ? resolve() : reject(new Error(stderr || `latexmk exited with code ${code}`)))
  })
}

export async function GET(request: Request, props: { params: Promise<{ fpl_id: string }> }) {
  const params = await props.params;
  return routeResponse(async () => {
    const query = new URL(request.url).searchParams
    const format = query.get('format') ?? 'zip'
    const metadataRaw = query.get('metadata')
    const metadata = metadataRaw ? JSON.parse(metadataRaw) : undefined
    const fpl = await fpprg.getFPL(params.fpl_id)
    if (!fpl) throw new NotFoundError()
    const session = await getRequestSessionWithId(request)
    const user = session?.user ? await getUser(session.user.id) : undefined
    const orcid = user ? await getUserOrcidAccount(user.id) : undefined
    const files = await FPL2TEX({
      krg: await getKRG(), fpl, metadata,
      author: user?.name ? {
        name: user.name,
        ...(user.email ? { email: user.email } : {}),
        ...(user.affiliation ? { affiliation: user.affiliation } : {}),
        ...(orcid ? { orcid: orcid.providerAccountId } : {}),
      } : undefined,
    })
    const names = Object.keys(files)
    if (format === 'zip') {
      const precomputed = names.find(name => name.endsWith('.zip'))
      const body = precomputed
        ? await files[precomputed] as ArrayBuffer
        : await (async () => {
          const zip = new JSZip()
          for (const [name, value] of Object.entries(files)) zip.file(name, value)
          return zip.generateAsync({ type: 'arraybuffer' })
        })()
      return new Response(body, { headers: {
        'Content-Type': 'application/zip',
        ...(precomputed ? { 'Content-Disposition': `attachment; filename="${precomputed}"` } : {}),
      } })
    }
    const precomputed = names.find(name => name.endsWith('.pdf'))
    if (precomputed) return new Response(await files[precomputed] as ArrayBuffer, { headers: {
      'Content-Type': 'application/pdf',
      'Content-Disposition': `attachment; filename="${precomputed}"`,
    } })
    const tmp = await fs.mkdtemp(path.join(os.tmpdir(), 'tex-'))
    try {
      for (const [name, pending] of Object.entries(files)) {
        const filename = path.join(tmp, name)
        await fs.mkdir(path.dirname(filename), { recursive: true })
        const value = await pending
        await fs.writeFile(filename, typeof value === 'string' ? value : new Uint8Array(value))
      }
      const main = names.find(name => name.endsWith('.tex') && !name.includes('/'))
      if (!main) return new Response('No main .tex file found in generated files', { status: 500 })
      await latexmk(path.join(tmp, main), tmp)
      const pdfName = main.replace(/\.tex$/, '.pdf')
      return new Response(await fs.readFile(path.join(tmp, pdfName)), { headers: {
        'Content-Type': 'application/pdf',
        'Content-Disposition': `attachment; filename="${pdfName}"`,
      } })
    } finally {
      await fs.rm(tmp, { recursive: true, force: true })
    }
  });
}
