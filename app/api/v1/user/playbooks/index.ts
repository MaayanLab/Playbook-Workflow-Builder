import { z } from 'zod'
import { API } from '@/spec/api'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { NotFoundError, UnauthorizedError } from '@/spec/error'
import db from '@/app/db'
import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import * as dict from '@/utils/dict'
import { IdOrPlaybookMetadataC, IdOrCellMetadataC } from '@/core/FPPRG'
import tsvector, { TSVector } from '@/utils/tsvector'

import publicPlaybooks from '@/app/public/playbooksDemo'
const playbook_tsvectors: Record<string, TSVector> = {}
publicPlaybooks.forEach(playbook => {
  playbook_tsvectors[playbook.id] = tsvector([
    playbook.label,
    playbook.description,
    ...playbook.inputs.flatMap(input => [
      input.meta.label,
      input.meta.description,
    ]),
    ...playbook.outputs.flatMap(output => [
      output.meta.label,
      output.meta.description,
    ]),
    ...playbook.dataSources,
  ].join(' '))
})

export const PublicPlaybooks = API('/api/v1/public/playbooks')
  .query(z.object({
    search: z.string().optional(),
    inputs: z.string().optional(),
    outputs: z.string().optional(),
    dataSources: z.string().optional(),
    skip: z.number().optional(),
    limit: z.number().optional(),
  }))
  .call(async (props, req, res) => {
    const search = props.query.search
    const inputs = props.query.inputs ? props.query.inputs.split(', ') : undefined
    const outputs = props.query.outputs ? props.query.outputs.split(', ') : undefined
    const dataSources = props.query.dataSources ? props.query.dataSources.split(', ') : undefined
    const skip = props.query.skip ?? 0
    const limit = props.query.limit ?? 50
    let playbooks = publicPlaybooks
    if (inputs) playbooks = playbooks.filter(playbook => !inputs.some(spec => !playbook.inputs.map(t=>t.spec as string).includes(spec)))
    if (outputs) playbooks = playbooks.filter(playbook => !outputs.some(spec => !playbook.outputs.map(t=>t.spec as string).includes(spec)))
    if (dataSources) playbooks = playbooks.filter(playbook => !dataSources.some(dataSource => !playbook.dataSources.includes(dataSource)))
    if (search) {
      const search_tsvector = tsvector(search)
      const search_scores: Record<string, number> = {}
      playbooks.forEach(playbook => {
        search_scores[playbook.id] = playbook_tsvectors[playbook.id]?.intersect(search_tsvector).size
      })
      playbooks = playbooks.filter(playbook => search_scores[playbook.id] > 0)
      playbooks.sort((a, b) => search_scores[b.id] - search_scores[a.id])
    }
    playbooks = playbooks.slice(skip, limit)
    const userPlaybooks = await db.objects.user_playbook.findMany({
      where: {
        playbook: {
          in: playbooks.map(({ id }) => id),
        },
      },
    })
    const userPlaybookLookup = dict.init(userPlaybooks.map((playbook) => ({ key: playbook.playbook, value: playbook })))
    const results = playbooks.map(({ inputs, outputs, dataSources, ...playbook }) => {
      return {
        ...playbook,
        inputs: inputs.map(t => t.spec as string).join(', '),
        outputs: outputs.map(t => t.spec as string).join(', '),
        dataSources: dataSources.join(', '),
        clicks: userPlaybookLookup[playbook.id]?.clicks || 0,
        disabled: userPlaybookLookup[playbook.id] === undefined,
      }
    })
    if (!search) results.sort((a, b) => b.clicks - a.clicks)
    return results
  })
  .build()

export const PublicUserPlaybooks = API('/api/v1/public/user/playbooks')
  .query(z.object({
    search: z.string().optional(),
    inputs: z.string().optional(),
    outputs: z.string().optional(),
    skip: z.number().optional(),
    limit: z.number().optional(),
  }))
  .call(async (props, req, res) => {
    const search = props.query.search
    const inputs = props.query.inputs ? props.query.inputs.split(', ') : undefined
    const outputs = props.query.outputs ? props.query.outputs.split(', ') : undefined
    const skip = props.query.skip ?? 0
    const limit = props.query.limit ?? 50
    // TODO: filter in DB
    let playbooks = await db.objects.user_playbook.findMany({
      where: {
        public: true,
      },
      orderBy: {
        clicks: 'desc',
      },
      skip,
      take: limit,
    })
    if (search) playbooks = playbooks.filter(playbook =>
      (playbook.title||'').includes(search)
      || (playbook.description||'').includes(search)
    )
    if (inputs) playbooks = playbooks.filter(playbook => !inputs.some(spec => !(playbook.inputs||'').split(', ').includes(spec)))
    if (outputs) playbooks = playbooks.filter(playbook => !outputs.some(spec => !(playbook.outputs||'').split(', ').includes(spec)))
    return playbooks
  })
  .build()

export const UserPlaybooks = API('/api/v1/user/playbooks')
  .query(z.object({
    skip: z.number().optional().transform(v => v ?? 0),
    limit: z.number().optional().transform(v => v ?? 50),
  }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const playbooks = await db.objects.user_playbook.findMany({
      where: {
        user: session.user.id,
      },
    })
    return playbooks
  })
  .build()

export const UserPlaybook = API('/api/v1/user/playbooks/[id]')
  .query(z.object({
    id: z.string(),
  }))
  .call(async (inputs, req, res) => {
    const fpl = await fpprg.getFPL(inputs.query.id)
    if (fpl === undefined) throw new NotFoundError()
    const metapath = fpl.resolve().map(fpl => fpl.toJSON())
    // TODO: move this logic into a db trigger
    // TODO: address playbook republishing
    const [publicUserPlaybook, ..._] = await db.objects.user_playbook.findMany({
      where: {
        playbook: inputs.query.id,
        public: true,
      },
    })
    if (publicUserPlaybook) {
      await db.objects.user_playbook.update({
        where: {
          id: publicUserPlaybook.id,
        },
        data: {
          clicks: publicUserPlaybook.clicks+1,
        }
      })
    }
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) return { metapath }
    const userPlaybook = await db.objects.user_playbook.findUnique({
      where: {
        user: session.user.id,
        playbook: inputs.query.id,
      },
    })
    if (!userPlaybook) return { metapath }
    return { metapath, userPlaybook: { public: userPlaybook.public } }
  })
  .build()

export const UpdateUserPlaybook = API('/api/v1/user/playbooks/[id]/update')
  .query(z.object({
    id: z.string(),
  }))
  .body(z.object({
    user_playbook: z.object({ public: z.boolean() }),
    playbook_metadata: IdOrPlaybookMetadataC,
    cell_metadata: z.record(z.string(), IdOrCellMetadataC),
  }))
  .call(async (inputs, req, res) => {
    let fpl = await fpprg.getFPL(inputs.query.id)
    if (fpl === undefined) throw new NotFoundError()
    fpl = fpl.rebaseCellMetadata(
      dict.init(
        await Promise.all(dict.items(inputs.body.cell_metadata)
          .map(async ({ key, value }) => ({ key, value: await fpprg.resolveCellMetadata(value) })))
      )
    )
    if (fpl === undefined) throw new NotFoundError()
    fpl = fpl.rebasePlaybookMetadata(fpl, await fpprg.resolvePlaybookMetadata(inputs.body.playbook_metadata)).head
    if (fpl === undefined) throw new NotFoundError()
    fpl = await fpprg.upsertFPL(fpl)
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const playbookInputs = fpl.resolve().map(cell => krg.getProcessNode(cell.process.type)).filter(node => 'prompt' in node).map(node => node.output.spec).join(', ')
    const playbookOutputs = fpl.resolve().filter(cell => cell.cell_metadata?.data_visible === true).map(cell => krg.getProcessNode(cell.process.type).output.spec).join(', ')
    await db.objects.user_playbook.delete({
      where: {
        user: session.user.id,
        playbook: inputs.query.id,
      }
    })
    await db.objects.user_playbook.upsert({
      where: {
        user: session.user.id,
        playbook: fpl.id,
      },
      create: {
        user: session.user.id,
        playbook: fpl.id,
        title: fpl.playbook_metadata?.title,
        description: fpl.playbook_metadata?.description,
        public: inputs.body.user_playbook.public,
        inputs: playbookInputs,
        outputs: playbookOutputs,
      },
      update: {
        title: fpl.playbook_metadata?.title,
        description: fpl.playbook_metadata?.description,
        public: inputs.body.user_playbook.public,
        inputs: playbookInputs,
        outputs: playbookOutputs,
      },
    })
    return fpl.id
  })
  .build()

export const PublishUserPlaybook = API('/api/v1/user/playbooks/[id]/publish')
  .query(z.object({
    id: z.string(),
  }))
  .body(z.object({
    public: z.boolean(),
  }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    return await db.objects.user_playbook.update({
      where: {
        user: session.user.id,
        playbook: inputs.query.id,
      },
      data: {
        public: inputs.body.public,
      },
    })
  })
  .build()

export const DeleteUserPlaybook = API('/api/v1/user/playbooks/[id]/delete')
  .query(z.object({
    id: z.string(),
  }))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const fpl = await fpprg.getFPL(inputs.query.id)
    if (fpl === undefined) throw new NotFoundError()
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const userPlaybook = await db.objects.user_playbook.delete({
      where: {
        user: session.user.id,
        playbook: inputs.query.id,
      },
    })
    return userPlaybook
  })
  .build()
