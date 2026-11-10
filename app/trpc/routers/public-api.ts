import { z } from 'zod'
import db from '@/app/db'
import fpprg from '@/app/fpprg'
import { getKRG } from '@/app/server-krg'
import packageJson from '@/package.json'
import { getUser, getUserOrcidAccount } from '@/app/extensions/auth/user'
import storage, { storageKey } from '@/components/core/file/server/storage'
import { drsUrl, resolveDrsObject } from '@/components/core/file/server/drs'
import FPL2BCO from '@/core/fpl2bco'
import fpl2json from '@/core/fpl2json'
import { cwl_for_playbook } from '@/core/fpl2cwl'
import { fpl2ro_crate_metadata } from '@/core/fpl2ro-crate'
import * as dict from '@/utils/dict'
import { NotFoundError } from '@/spec/error'
import { DataC, IdOrPlaybookMetadataC, IdOrProcessC } from '@/core/FPPRG'
import { publicProcedure, router } from '../init'

const JsonValue: z.ZodType<unknown> = z.lazy(() => z.union([
  z.string(), z.number(), z.boolean(), z.null(), z.array(JsonValue), z.record(z.string(), JsonValue),
]))

/** Match the JSON representation returned on the wire (and remove undefined object fields). */
function asJson(value: unknown): z.infer<typeof JsonValue> {
  return JSON.parse(JSON.stringify(value)) as z.infer<typeof JsonValue>
}

const FplIdInput = z.object({ fpl_id: z.string() })

const DrsServiceInfo = z.object({
  id: z.string(),
  name: z.string(),
  type: z.object({ group: z.string(), artifact: z.string(), version: z.string() }),
  description: z.string(),
  organization: z.object({ name: z.string(), url: z.string() }),
  contactUrl: z.string(),
  documentationUrl: z.string(),
  createdAt: z.string(),
  updatedAt: z.string(),
  environment: z.string(),
  version: z.string(),
})

const DrsObject = z.object({
  id: z.string(),
  name: z.string(),
  self_uri: z.string(),
  size: z.number(),
  created_time: z.union([z.string(), z.date()]),
  checksums: z.array(z.object({ type: z.string(), checksum: z.string() })),
  access_methods: z.array(z.union([
    z.object({ type: z.string(), access_id: z.string() }),
    z.object({ type: z.string(), access_url: z.string() }),
  ])),
})

export async function getWorkflowExport(fplId: string) {
  const fpl = await fpprg.getFPL(fplId)
  if (!fpl) throw new NotFoundError()
  const data: Record<string, unknown> = {}
  const workflow = []
  for (const el of fpl.resolve()) {
    workflow.push({
      id: el.process.id,
      type: el.process.type,
      inputs: !dict.isEmpty(el.process.inputs)
        ? dict.init(dict.items(el.process.inputs).map(({ key, value }) => ({ key, value: { id: value.id } })))
        : undefined,
      data: el.process.data ? { id: el.process.data.id } : undefined,
    })
    if (el.process.data) {
      const { id: _, ...dataNode } = el.process.data.toJSON()
      data[el.process.data.id] = dataNode
    }
  }
  return { data, workflow, metadata: fpl.playbook_metadata?.toJSON() }
}

/**
 * Where an upload can be read from: directly from storage when it can give out (time limited)
 *  urls, e.g. presigned S3 urls, or from its source if it's remote, otherwise from this server.
 */
async function uploadAccessUrl(upload: { id: string, url: string, filename: string }) {
  if (upload.url.startsWith('storage://')) {
    const backend = storage()
    if (backend.presignedGet) return await backend.presignedGet(storageKey(upload.url), { filename: upload.filename })
  } else if (upload.url.startsWith('https://')) {
    return upload.url
  }
  return `${process.env.PUBLIC_URL || ''}/api/ga4gh/drs/v1/objects/${upload.id}/access/https/data`
}

export const publicApiRouter = router({
  drsServiceInfo: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/ga4gh/drs/v1/service-info', tags: ['GA4GH DRS'],
      summary: 'Get DRS service information',
    } })
    .input(z.object({}))
    .output(DrsServiceInfo)
    .query(() => ({
      id: 'cloud.playbook-workflow-builder',
      name: 'Playbook Workflow Builder',
      type: { group: 'org.ga4gh', artifact: 'drs', version: '1.0.0' },
      description: 'This service provides DRS capabilities for the playbook workflow builder.',
      organization: { name: 'Playbook Partnership', url: 'https://playbook-workflow-builder.cloud' },
      contactUrl: 'mailto:avi.maayan@mssm.edu',
      documentationUrl: 'https://github.com/MaayanLab/Playbook-Workflow-Builder',
      createdAt: '2023-04-20T12:20:15Z',
      updatedAt: '2023-09-05T12:20:15Z',
      environment: 'prod',
      version: packageJson.version,
    })),

  drsObject: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/ga4gh/drs/v1/objects/{object_id}', tags: ['GA4GH DRS'],
      summary: 'Get a DRS object', errorResponses: [404],
    } })
    .input(z.object({ object_id: z.string() }))
    .output(DrsObject)
    .query(async ({ input }) => {
      const upload = await resolveDrsObject(input.object_id)
      if (!upload) throw new NotFoundError()
      return {
        id: upload.id,
        name: upload.filename,
        self_uri: drsUrl(upload.id),
        size: upload.size,
        created_time: upload.created,
        checksums: [{ type: 'sha-256', checksum: upload.sha256 }],
        access_methods: [
          { type: 'https', access_id: 'https' },
          { type: 'https', access_url: `${process.env.PUBLIC_URL || ''}/api/ga4gh/drs/v1/objects/${input.object_id}/access/https/data` },
        ],
      }
    }),

  drsAccess: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/ga4gh/drs/v1/objects/{object_id}/access/{access_id}', tags: ['GA4GH DRS'],
      summary: 'Get a DRS access URL', errorResponses: [404],
    } })
    .input(z.object({ object_id: z.string(), access_id: z.string() }))
    .output(z.object({ url: z.string() }))
    .query(async ({ input }) => {
      if (input.access_id !== 'https') throw new NotFoundError()
      const upload = await resolveDrsObject(input.object_id)
      if (!upload) throw new NotFoundError()
      return { url: await uploadAccessUrl(upload) }
    }),

  bco: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/bco/{fpl_id}', tags: ['Workflows'], summary: 'Export a playbook as a BCO',
      errorResponses: [404],
    } })
    .input(FplIdInput.extend({ metadata: z.string().optional() }))
    .output(JsonValue)
    .query(async ({ input, ctx }) => {
      const fpl = await fpprg.getFPL(input.fpl_id)
      if (!fpl) throw new NotFoundError()
      const metadata = input.metadata
        ? z.object({ title: z.string().optional(), description: z.string().optional() }).parse(JSON.parse(input.metadata))
        : undefined
      const session = await ctx.getSession()
      const user = session?.user ? await getUser(session.user.id) : undefined
      const userOrcidAccount = user ? await getUserOrcidAccount(user.id) : undefined
      const krg = await getKRG()
      return asJson(await FPL2BCO({
        krg, fpl, metadata,
        author: user?.name ? {
          name: user.name,
          ...dict.filter({
            email: user.email,
            affiliation: user.affiliation ?? undefined,
            orcid: userOrcidAccount?.providerAccountId,
          }),
        } : undefined,
      }))
    }),

  workflowCreate: publicProcedure
    .meta({ openapi: {
      method: 'POST', path: '/db/fpl', tags: ['Workflows'],
      summary: 'Create a playbook from a workflow specification, returning its id',
      description: 'The body uses the same structure as `/db/fpl/{fpl_id}/export`.',
      errorResponses: [400],
    } })
    .input(z.object({
      data: z.record(z.string(), DataC).optional(),
      workflow: z.array(IdOrProcessC),
      metadata: IdOrPlaybookMetadataC.optional(),
    }))
    .output(z.string())
    .mutation(async ({ input }) => (await fpprg.resolveFPL(input)).id),

  workflow: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/db/fpl/{fpl_id}', tags: ['Workflows'],
      summary: 'Get the steps of a playbook as registered in the database', errorResponses: [404],
    } })
    .input(FplIdInput)
    .output(JsonValue)
    .query(async ({ input }) => {
      const fpl = await fpprg.getFPL(input.fpl_id)
      if (!fpl) throw new NotFoundError()
      return asJson(fpl.resolve().map(step => step.toJSON()))
    }),

  workflowOutput: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/db/fpl/{fpl_id}/output', tags: ['Workflows'],
      summary: 'Get the steps of a playbook along with their outputs', errorResponses: [404, 504],
    } })
    .input(FplIdInput)
    .output(JsonValue)
    .query(async ({ input }) => {
      const fpl = await fpprg.getFPL(input.fpl_id)
      if (!fpl) throw new NotFoundError()
      return asJson(await Promise.all(fpl.resolve().map(step => step.toJSONWithOutput())))
    }),

  process: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/db/process/{process_id}', tags: ['Workflows'],
      summary: 'Get a process as registered in the database', errorResponses: [404],
    } })
    .input(z.object({ process_id: z.string() }))
    .output(JsonValue)
    .query(async ({ input }) => {
      const process = await fpprg.getProcess(input.process_id)
      if (!process) throw new NotFoundError()
      return asJson(process.toJSON())
    }),

  processOutput: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/db/process/{process_id}/output', tags: ['Workflows'],
      summary: 'Get the output of a process, resolving it if necessary',
      description: 'Responds with 504 if the output is not ready yet, retry until it is.',
      errorResponses: [404, 504],
    } })
    .input(z.object({ process_id: z.string() }))
    .output(JsonValue)
    .query(async ({ input }) => {
      const process = await fpprg.getProcess(input.process_id)
      if (!process) throw new NotFoundError()
      const output = await process.output()
      return asJson(output ? output.toJSON() : null)
    }),

  workflowExport: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/db/fpl/{fpl_id}/export', tags: ['Workflows'],
      summary: 'Export a playbook workflow', errorResponses: [404],
    } })
    .input(FplIdInput)
    .output(JsonValue)
    .query(async ({ input }) => asJson(await getWorkflowExport(input.fpl_id))),

  outputExport: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/db/fpl/{fpl_id}/output/export', tags: ['Workflows'],
      summary: 'Export a playbook with outputs', errorResponses: [404],
    } })
    .input(FplIdInput)
    .output(JsonValue)
    .query(async ({ input }) => {
      const fpl = await fpprg.getFPL(input.fpl_id)
      if (!fpl) throw new NotFoundError()
      return asJson(await fpl2json({ fpl }))
    }),

  cwl: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/v1/cwl/{fpl_id}', tags: ['Workflows'],
      summary: 'Export CWL as JSON or list generated files', errorResponses: [404],
    } })
    .input(FplIdInput.extend({ format: z.enum(['json', 'list']) }))
    .output(JsonValue)
    .query(async ({ input }) => {
      const fpl = await fpprg.getFPL(input.fpl_id)
      if (!fpl) throw new NotFoundError()
      const krg = await getKRG()
      const files = await cwl_for_playbook({ krg, fpl })
      return asJson(input.format === 'list' ? dict.keys(files) : files)
    }),

  roCrate: publicProcedure
    .meta({ openapi: {
      method: 'GET', path: '/v1/ro-crate/{fpl_id}', tags: ['Workflows'],
      summary: 'Export RO-Crate metadata as JSON', errorResponses: [404],
    } })
    .input(FplIdInput.extend({ format: z.literal('json') }))
    .output(JsonValue)
    .query(async ({ input }) => {
      const fpl = await fpprg.getFPL(input.fpl_id)
      if (!fpl) throw new NotFoundError()
      const krg = await getKRG()
      return asJson(await fpl2ro_crate_metadata({ krg, fpl }))
    }),
})
