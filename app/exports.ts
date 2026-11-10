/**
 * Playbook exports (LaTeX bundles, report PDFs, RO-Crates, CWL bundles) can take longer to produce than
 *  the ingress allows a request to take. So rather than producing them in the request,
 *  they're submitted as jobs which are processed by the same workers that resolve processes,
 *  and their progress is reported over SSE (`exports.watch`). Since they produce files which
 *  are kept, exports require signing in: an export and its file (one of the user's uploads,
 *  listed on their account page) belong to the user who requested it.
 *
 * The synchronous routes (`/api/v1/tex`, `/api/v1/ro-crate`, `/api/v1/cwl`) remain for API users and produce
 *  the same files with `generateExport`.
 */
import { spawn } from 'child_process'
import fs from 'fs/promises'
import os from 'os'
import path from 'path'
import JSZip from 'jszip'
import YAML from 'yaml'
import { z } from 'zod'
import db, { bus } from '@/app/db'
import { uuid } from '@/db'
import { selectUserUploadComplete } from '@/db/queries'
import type { FplExportStatus } from '@/db/schema'
import type FPPRG from '@/core/FPPRG'
import type { FPL } from '@/core/FPPRG'
import type KRG from '@/core/KRG'
import FPL2TEX from '@/core/fpl2tex'
import { fpl2ro_crate } from '@/core/fpl2ro-crate'
import { cwl_for_playbook } from '@/core/fpl2cwl'
import type { Author } from '@/core/common'
import { getUser, getUserOrcidAccount } from '@/app/extensions/auth/user'
import { uploadFile } from '@/components/core/file/server/upload'
import { NotFoundError } from '@/spec/error'

export const ExportFormatC = z.enum(['tex-zip', 'tex-pdf', 'ro-crate-zip', 'cwl-zip'])
export type ExportFormat = z.infer<typeof ExportFormatC>

export const ExportMetadataC = z.object({
  title: z.string().nullish().transform(v => v ?? undefined),
  description: z.string().nullish().transform(v => v ?? undefined),
})

const AuthorC = z.object({
  name: z.string(),
  affiliation: z.string().optional(),
  email: z.string().optional(),
  orcid: z.string().optional(),
})

const ExportParamsC = z.object({
  metadata: ExportMetadataC.optional(),
  author: AuthorC.optional(),
})
export type ExportParams = z.infer<typeof ExportParamsC>

export type ExportFile = { filename: string, contentType: string, data: Uint8Array<ArrayBuffer> }

/** The queue export jobs are submitted to */
const EXPORT_QUEUE = 'export-queue'
const exportEvent = (id: string) => `distributed:export:${id}`

/** The author of the report: the signed in user, with their ORCID if they've linked it */
export async function exportAuthor(userId: string | undefined): Promise<Author | undefined> {
  const user = userId ? await getUser(userId) : undefined
  if (!user?.name) return undefined
  const orcid = await getUserOrcidAccount(user.id)
  return {
    name: user.name,
    ...(user.email ? { email: user.email } : {}),
    ...(user.affiliation ? { affiliation: user.affiliation } : {}),
    ...(orcid ? { orcid: orcid.providerAccountId } : {}),
  }
}

function latexmk(main: string, cwd: string) {
  return new Promise<void>((resolve, reject) => {
    const proc = spawn('latexmk', ['-pdf', '-cd', '-lualatex', '-interaction=nonstopmode', '-f', main], { cwd })
    let stderr = ''
    proc.stdout.on('data', data => process.stdout.write(`[latexmk] ${data}`))
    proc.stderr.on('data', data => { stderr += data; process.stderr.write(`[latexmk] ${data}`) })
    proc.on('error', reject)
    proc.on('close', code => code === 0 ? resolve() : reject(new Error(stderr || `latexmk exited with code ${code}`)))
  })
}

const asBytes = (value: string | ArrayBuffer | Uint8Array): Uint8Array<ArrayBuffer> =>
  typeof value === 'string' ? new TextEncoder().encode(value)
  : new Uint8Array(value)

async function zipFiles(files: Record<string, Promise<string | ArrayBuffer> | string | ArrayBuffer>) {
  const zip = new JSZip()
  for (const [name, value] of Object.entries(files)) zip.file(name, await value)
  return new Uint8Array(await zip.generateAsync({ type: 'arraybuffer' }))
}

/**
 * Produce an export of `fpl`
 */
export async function generateExport({ krg, fpl, format, params }: { krg: KRG, fpl: FPL, format: ExportFormat, params: ExportParams }): Promise<ExportFile> {
  if (format === 'cwl-zip') {
    const files = await cwl_for_playbook({ krg, fpl })
    const zip = new JSZip()
    for (const [name, value] of Object.entries(files)) zip.file(name, YAML.stringify(value))
    return { filename: `cwl-${fpl.id}.zip`, contentType: 'application/zip', data: new Uint8Array(await zip.generateAsync({ type: 'arraybuffer' })) }
  }
  if (format === 'ro-crate-zip') {
    const files = await fpl2ro_crate({ krg, fpl })
    return { filename: `ro-crate-${fpl.id}.zip`, contentType: 'application/zip', data: await zipFiles(files) }
  }
  const files = await FPL2TEX({ krg, fpl, metadata: params.metadata, author: params.author })
  const names = Object.keys(files)
  // a report node produces the bundle & pdf itself, they're at the top level
  //  (unlike e.g. figures/1.pdf)
  const topLevel = (ext: string) => names.find(name => name.endsWith(ext) && !name.includes('/'))
  if (format === 'tex-zip') {
    const precomputed = topLevel('.zip')
    return {
      filename: precomputed ?? `${fpl.id}.zip`,
      contentType: 'application/zip',
      data: precomputed ? asBytes(await files[precomputed]) : await zipFiles(files),
    }
  }
  const precomputed = topLevel('.pdf')
  if (precomputed) return { filename: precomputed, contentType: 'application/pdf', data: asBytes(await files[precomputed]) }
  const main = topLevel('.tex')
  if (!main) throw new Error('No main .tex file found in generated files')
  const tmp = await fs.mkdtemp(path.join(os.tmpdir(), 'tex-'))
  try {
    for (const [name, pending] of Object.entries(files)) {
      const filename = path.join(tmp, name)
      await fs.mkdir(path.dirname(filename), { recursive: true })
      await fs.writeFile(filename, asBytes(await pending))
    }
    await latexmk(path.join(tmp, main), tmp)
    const filename = main.replace(/\.tex$/, '.pdf')
    return { filename, contentType: 'application/pdf', data: new Uint8Array(await fs.readFile(path.join(tmp, filename))) }
  } finally {
    await fs.rm(tmp, { recursive: true, force: true })
  }
}

export type ExportState = {
  id: string
  format: ExportFormat
  status: FplExportStatus
  error: string | null
  /** available once complete */
  download: { id: string, filename: string, size: number } | null
}

/** The state of export `id`, if it belongs to `user` */
export async function getExport(id: string, user: string): Promise<ExportState | undefined> {
  const row = await db.selectFrom('fpl_export')
    .select(['id', 'format', 'status', 'error', 'upload'])
    .where('id', '=', id)
    .where('user', '=', user)
    .executeTakeFirst()
  if (!row) return undefined
  // (the user may since have deleted the file)
  const upload = row.upload ? await selectUserUploadComplete(db)
    .where('user_upload.id', '=', row.upload)
    .where('user_upload.user', '=', user)
    .executeTakeFirst() : undefined
  return {
    id: row.id,
    format: row.format as ExportFormat,
    status: row.status,
    error: row.error,
    download: upload ? { id: upload.id, filename: upload.filename, size: upload.size } : null,
  }
}

/**
 * Queue an export of `fpl_id`, returning its id (watch it with `watchExport`)
 */
export async function requestExport({ fpl_id, format, params, user }: { fpl_id: string, format: ExportFormat, params: ExportParams, user: string }) {
  const id = uuid()
  await db.insertInto('fpl_export')
    .values({ id, fpl: fpl_id, user, format, params: JSON.stringify(params), status: 'queued' })
    .execute()
  await bus.send(EXPORT_QUEUE, { id })
  return id
}

async function updateExport(id: string, values: { status: FplExportStatus, error?: string | null, upload?: string | null }) {
  await db.updateTable('fpl_export')
    .set({ ...values, updated: new Date().toISOString() })
    .where('id', '=', id)
    .execute()
  bus.notify(exportEvent(id), { status: values.status })
}

/**
 * Report the state of export `id` each time it changes, until it's complete (or failed).
 */
export async function* watchExport(id: string, user: string, signal?: AbortSignal): AsyncGenerator<ExportState> {
  // listen before reading the state so that no change can be missed in between
  let dirty = true
  let wake: (() => void) | undefined
  const unlisten = bus.listen(exportEvent(id), () => { dirty = true; wake?.() })
  const onAbort = () => wake?.()
  signal?.addEventListener('abort', onAbort, { once: true })
  try {
    let last: string | undefined
    while (!signal?.aborted) {
      if (dirty) {
        dirty = false
        const state = await getExport(id, user)
        if (!state) throw new NotFoundError()
        const serialized = JSON.stringify(state)
        if (serialized !== last) {
          last = serialized
          yield state
        }
        if (state.status === 'complete' || state.status === 'error') return
      }
      // re-check periodically in case a notification was lost (e.g. a database reconnect)
      await new Promise<void>(resolve => {
        const timeout = setTimeout(() => { dirty = true; resolve() }, 30_000)
        wake = () => { clearTimeout(timeout); resolve() }
      })
      wake = undefined
    }
  } finally {
    unlisten()
    signal?.removeEventListener('abort', onAbort)
  }
}

/**
 * Process export jobs, run alongside the process resolution workers.
 *  Failures are recorded on the export rather than thrown (which would have the job retried).
 */
export function start_export_workers(getKRG: () => Promise<KRG>, fpprg: FPPRG, concurrency: number) {
  return bus.work(EXPORT_QUEUE, { concurrency }, async (job) => {
    const id = job.data.id
    const row = await db.selectFrom('fpl_export').selectAll().where('id', '=', id).executeTakeFirst()
    if (!row) {
      console.error(`export '${id}' not found`)
      return
    }
    if (row.status === 'complete') return
    console.debug(`exporting ${row.format} of ${row.fpl} (${id})..`)
    await updateExport(id, { status: 'running', error: null })
    try {
      const fpl = await fpprg.getFPL(row.fpl)
      if (!fpl) throw new NotFoundError()
      const file = await generateExport({
        krg: await getKRG(),
        fpl,
        format: ExportFormatC.parse(row.format),
        params: ExportParamsC.parse(JSON.parse(row.params)),
      })
      const tmp = await fs.mkdtemp(path.join(os.tmpdir(), 'pwb-export-'))
      try {
        const tmpFile = path.join(tmp, path.basename(file.filename))
        await fs.writeFile(tmpFile, file.data)
        // the file is one of the requester's uploads
        const upload = await uploadFile({
          url: `file://${tmpFile.split(path.sep).join(path.posix.sep)}`,
          size: file.data.length,
          filename: file.filename,
        }, { user: { id: row.user } })
        await updateExport(id, { status: 'complete', upload: upload.id })
      } finally {
        await fs.rm(tmp, { recursive: true, force: true })
      }
      console.debug(`completed export ${id}`)
    } catch (error) {
      console.error(`export ${id} failed`, error)
      const message = error instanceof Error ? error.message : String(error)
      await updateExport(id, { status: 'error', error: message.slice(0, 4000) })
    }
  })
}
