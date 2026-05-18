import { z } from 'zod'
import { API } from '@/spec/api'
import FPL2TEX from '@/core/fpl2tex'
import fpprg from '@/app/fpprg'
import krg from '@/app/krg'
import { NotFoundError } from '@/spec/error'
import JSZip from 'jszip'
import * as dict from '@/utils/dict'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import db from '@/app/db'
import { exec,spawn } from 'child_process'
import fs from 'fs'
import os from 'os'
import path from 'path'

export const TeXForPlaybook = API.get('/api/v1/tex/[fpl_id]')
  .query(z.object({
    fpl_id: z.string(),
    format: z.enum(['zip','pdf']).default('zip'),
    metadata: z.object({
      title: z.string().optional(),
      description: z.string().optional(),
    }).optional(),
  }))
  .call(async (inputs, req, res) => {
    const fpl_id = inputs.query.fpl_id
    const metadata = inputs.query.metadata
    
    const fpl = await fpprg.getFPL(fpl_id)
    if (fpl === undefined) throw new NotFoundError()
    const session = await getServerSessionWithId(req, res)
    const user = (session && session.user) ? (await db.objects.user.findUnique({ where: { id: session.user.id } })) : undefined
    // @ts-ignore
    const userOrcidAccount = user ? await db.objects.account.findUnique({ where: { userId: user.id, provider: 'orcid' } }) : undefined
    const tex_files = await FPL2TEX({
      krg,
      fpl,
      metadata,
      author: user && user.name !== undefined ? {
        name: user.name,
        ...dict.filter({
          email: user.email,
          affiliation: user.affiliation,
          orcid: userOrcidAccount ? userOrcidAccount.providerAccountId : undefined,
        }),
      } : undefined,
    })
    if (inputs.query.format === 'zip') {
      const zip = new JSZip()
      dict.items(tex_files).forEach(({ key, value }) => {
        zip.file(key, value)
      })
      zip.generateNodeStream({
        type: 'nodebuffer',
        streamFiles: true,
      }).pipe(res)
    } else if (inputs.query.format === 'pdf') {
      const tmpDir = fs.mkdtempSync(path.join(os.tmpdir(), 'tex-'))
      try {
        for (const { key, value } of dict.items(tex_files)) {
          const filePath = path.join(tmpDir, key)
          fs.mkdirSync(path.dirname(filePath), { recursive: true })
          const fileValue = await value
          fs.writeFileSync(filePath,typeof fileValue === 'string' ? fileValue : Buffer.from(fileValue))
        
        }
        const mainTex = dict.keys(tex_files).find(k => k.endsWith('.tex') && !k.includes('/'))
        if (!mainTex) {
          res.status(500).send('No main .tex file found in generated files')
          return
        }
        const mainTexPath = path.join(tmpDir, mainTex)
        const files = fs.readdirSync(tmpDir)
        await new Promise<void>((resolve, reject) => {
          const proc = spawn(
            'latexmk',
            ['-pdf', '-cd', '-lualatex', '-interaction=nonstopmode', '-halt-on-error', mainTexPath],
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
    
        const pdfPath = mainTexPath.replace(/\.tex$/, '.pdf')
        const pdfBuffer = fs.readFileSync(pdfPath)
        res.setHeader('Content-Type', 'application/pdf')
        res.setHeader('Content-Disposition', `attachment; filename="${mainTex.replace(/\.tex$/, '.pdf')}"`)
        res.send(pdfBuffer)
      } finally {
        fs.rmSync(tmpDir, { recursive: true, force: true })
      }
    }
  })
  .build()
