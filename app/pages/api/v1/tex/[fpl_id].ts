import { z } from 'zod'
import FPL2TEX from '@/core/fpl2tex'
import fpprg from '@/app/fpprg'
import { getKRG } from '@/app/server-krg'
import { NotFoundError, UnsupportedMethodError } from '@/spec/error'
import JSZip from 'jszip'
import * as dict from '@/utils/dict'
import handler from '@/utils/next-rest'
import { getServerSessionWithId } from '@/app/extensions/auth/server'
import { getUser, getUserOrcidAccount } from '@/app/extensions/auth/user'
import { spawn } from 'child_process'
import fs from 'fs/promises'
import os from 'os'
import path from 'path'

function tryJsonParse(value: unknown) {
  if (typeof value !== 'string') return value
  try { return JSON.parse(value) } catch { return value }
}

const QueryType = z.object({
  fpl_id: z.string(),
  format: z.enum(['zip','pdf']).default('zip'),
  metadata: z.preprocess(tryJsonParse, z.object({
    title: z.string().optional(),
    description: z.string().optional(),
  }).optional()),
})

export default handler(async (req, res) => {
  if (req.method !== 'GET') throw new UnsupportedMethodError(req.method)
  const { fpl_id, format, metadata } = QueryType.parse(req.query)
  
  const fpl = await fpprg.getFPL(fpl_id)
  if (fpl === undefined) throw new NotFoundError()
  const session = await getServerSessionWithId(req, res)
  const user = (session && session.user) ? (await getUser(session.user.id)) : undefined
  const userOrcidAccount = user ? await getUserOrcidAccount(user.id) : undefined
  const krg = await getKRG()
  const tex_files = await FPL2TEX({
    krg,
    fpl,
    metadata,
    author: user && user.name ? {
      name: user.name,
      ...dict.filter({
        email: user.email,
        affiliation: user.affiliation ?? undefined,
        orcid: userOrcidAccount ? userOrcidAccount.providerAccountId : undefined,
      }),
    } : undefined,
  })
  if (format === 'zip') {
    const precomputedZip = dict.keys(tex_files).find(k => k.endsWith('.zip'))
    if (precomputedZip) {
      res.setHeader('Content-Type', 'application/zip')
      res.setHeader('Content-Disposition', `attachment; filename="${precomputedZip}"`)
      const buffer = await tex_files[precomputedZip]
      res.end(Buffer.from(buffer as ArrayBuffer))
      return
    }
    const zip = new JSZip()
    dict.items(tex_files).forEach(({ key, value }) => {
      zip.file(key, value)
    })
    zip.generateNodeStream({
      type: 'nodebuffer',
      streamFiles: true,
    }).pipe(res)
  } else if (format === 'pdf') {
    const precomputedPdf = dict.keys(tex_files).find(k => k.endsWith('.pdf'))
    if (precomputedPdf) {
      res.setHeader('Content-Type', 'application/pdf')
      res.setHeader('Content-Disposition', `attachment; filename="${precomputedPdf}"`)
      const buffer = await tex_files[precomputedPdf]
      res.end(Buffer.from(buffer as ArrayBuffer))
      return
    }
    const tmpDir = await fs.mkdtemp(path.join(os.tmpdir(), 'tex-'))
    try {
      res.setHeader('Content-Type', 'application/pdf')
      await Promise.allSettled(dict.items(tex_files).map(async ({ key, value }) => {
        const filePath = path.join(tmpDir, key)
        await fs.mkdir(path.dirname(filePath), { recursive: true })
        const fileValue = await value
        await fs.writeFile(filePath,typeof fileValue === 'string' ? fileValue : new DataView(fileValue))
      }))
      const mainTex = dict.keys(tex_files).find(k => k.endsWith('.tex') && !k.includes('/'))
      if (!mainTex) {
        res.status(500).send('No main .tex file found in generated files')
        return
      }
      const mainTexPath = path.join(tmpDir, mainTex)
      const pdfPath = mainTexPath.replace(/\.tex$/, '.pdf')
      res.setHeader('Content-Disposition', `attachment; filename="${mainTex.replace(/\.tex$/, '.pdf')}"`)
      try {
        await new Promise<void>((resolve, reject) => {
          const proc = spawn(
            'latexmk',
            ['-pdf', '-cd', '-lualatex', '-interaction=nonstopmode', '-f', mainTexPath],
            { cwd: tmpDir }
          )
        
          let stderr = ''
          proc.stdout.on('data', (data) => process.stdout.write(`[latexmk] ${data}`))
          proc.stderr.on('data', (data) => {
            process.stderr.write(`[latexmk] ${data}`)
            stderr += data
          })
        
          proc.on('close', (code) => {
            if (code !== 0) reject(new Error(stderr || `latexmk exited with code ${code}`))
            else resolve()
          })
        })
      } catch (e) {
        console.error(e)
      } finally {
        const pdfBuffer = await fs.open(pdfPath)
        pdfBuffer.createReadStream({ autoClose: true }).pipe(res)
      }
    } catch (e) {
      console.error(e)
      res.status(500)
      res.end()
    } finally {
      await fs.rm(tmpDir, { recursive: true, force: true })
    }
  }
})
