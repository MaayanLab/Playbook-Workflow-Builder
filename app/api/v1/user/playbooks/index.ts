import { z } from 'zod'
import { API } from '@/spec/api'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { NotFoundError, UnauthorizedError } from '@/spec/error'
import db from '@/app/db'
import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import * as dict from '@/utils/dict'
import { IdOrPlaybookMetadataC, IdOrCellMetadataC } from '@/core/FPPRG'

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
    const limit = props.query.limit ?? 10
    // TODO: filter in DB
    let playbooks = await db.objects.user_playbook.findMany({
      where: {
        public: true,
      },
    })
    if (search) playbooks = playbooks.filter(playbook =>
      (playbook.title||'').includes(search)
      || (playbook.description||'').includes(search)
    )
    if (inputs) playbooks = playbooks.filter(playbook => !inputs.some(spec => !(playbook.inputs||'').split(', ').includes(spec)))
    if (outputs) playbooks = playbooks.filter(playbook => !outputs.some(spec => !(playbook.outputs||'').split(', ').includes(spec)))
    return playbooks.slice(skip, skip + limit)
  })
  .build()

export const UserPlaybooks = API('/api/v1/user/playbooks')
  .query(z.object({
    skip: z.number().optional().transform(v => v ?? 0),
    limit: z.number().optional().transform(v => v ?? 10),
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
