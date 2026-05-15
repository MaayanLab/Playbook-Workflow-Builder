import { z } from 'zod'
import { API } from '@/spec/api'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { NotFoundError, UnauthorizedError } from '@/spec/error'
import db from '@/app/db'
import * as schema from "@/db"
import krg from '@/app/krg'
import fpprg from '@/app/fpprg'
import * as dict from '@/utils/dict'
import { IdOrPlaybookMetadataC, IdOrCellMetadataC } from '@/core/FPPRG'
import { tsvector, tsvector_intersect } from '@/utils/tsvector'
import cache from '@/utils/global_cache'
import { embed, indexSoon } from '@/core/semantic-search'

import publicPlaybooks from '@/app/public/playbooksDemo'
import { fpl_expand } from '@/core/common'
import { PgDatabase } from '@/utils/orm/pg'
import openai from '@/app/extensions/openai'
const with_public_playbooks = cache('with_public_playbooks', async () => {
  await Promise.all(publicPlaybooks.map(async (playbook) => {
    try {
      const fpl = await fpprg.upsertFPL(await fpprg.resolveFPL(playbook.workflow as any))
      const { story } = await fpl_expand({ krg, fpl })
      const storyText = story.ast.flatMap(part => part.tags.includes('abstract') ? [part.type === 'bibitem' ? '\n' : '', part.text] : []).join('')
      const playbookToInsert = {
        id: playbook.id,
        label: playbook.label,
        authors: playbook.authors.join('\t'),
        data_sources: playbook.dataSources.join('\t'),
        story: storyText,
        description: playbook.description,
        gpt_summary: playbook.gpt_summary,
        inputs: playbook.inputs.map(inp => inp.spec).join('\t'),
        outputs: playbook.outputs.map(inp => inp.spec).join('\t'),
        license: playbook.license,
        licenseUrl: playbook.licenseUrl,
        playbook: fpl.id,
        published: new Date(playbook.published),
        version: playbook.version,
      }
      await db.objects.published_playbook.upsert({
        where: { id: playbook.id },
        create: playbookToInsert,
        update: playbookToInsert,
      })
      await indexSoon({ db, fpl_id: fpl.id })
    } catch (e) {
      console.error(e)
    }
  }))
})

export const PublicPlaybooks = API.get('/api/v1/public/playbooks')
  .query(z.object({
    search: z.string().optional(),
    inputs: z.string().optional(),
    outputs: z.string().optional(),
    dataSources: z.string().optional(),
    skip: z.number().optional(),
    limit: z.number().optional(),
  }))
  .call(async (props) => {
    await with_public_playbooks
    const search = props.query.search
    const inputs = props.query.inputs ? props.query.inputs.split('\t') : undefined
    const outputs = props.query.outputs ? props.query.outputs.split('\t') : undefined
    const dataSources = props.query.dataSources ? props.query.dataSources.split('\t') : undefined
    const skip = props.query.skip ?? 0
    const limit = props.query.limit ?? 50
    if ('db' in db && db.db instanceof PgDatabase) {
      const embedding = search ? await embed({ openai: await openai, content: search }) : undefined
      const playbooks = (await db.db.raw(subst => `
        select "published_playbook".*
        from "published_playbook"
        ${search ? `left join fpl_embedding on "published_playbook".playbook = fpl_embedding.id` : ''}
        ${inputs || outputs || dataSources ? `where ${[
          inputs && `"published_playbook".inputs like any(${subst(inputs.map(s => `%${s}%`))})`,
          outputs && `"published_playbook".outputs like any(${subst(outputs.map(s => `%${s}%`))})`,
          dataSources && `"published_playbook".data_sources like any(${subst(dataSources.map(s => `%${s}%`))})`,
        ].filter(v => !!v).join(' and ')}` : ''}
        ${search ? `order by "fpl_embedding".embedding <-> ${subst(embedding)}` : `order by "clicks" desc`}
        offset ${subst(skip)}
        limit ${subst(limit)}
      `)).rows.map(row => schema.published_playbook.codec.decode(row))
      return playbooks
    } else {
      let playbooks = await db.objects.published_playbook.findMany({
        orderBy: {
          clicks: 'desc',
        },
        skip,
        take: limit,
      })
      if (inputs) playbooks = playbooks.filter(playbook => !inputs.some(spec => !playbook.inputs.split('\t').includes(spec)))
      if (outputs) playbooks = playbooks.filter(playbook => !outputs.some(spec => !playbook.outputs.split('\t').includes(spec)))
      if (dataSources) playbooks = playbooks.filter(playbook => !dataSources.some(spec => !playbook.data_sources?.split('\t').includes(spec)))
      if (search) {
        const search_tsvector = tsvector(search)
        const search_scores = dict.init(playbooks.map(playbook => {
          const playbook_tsvector = tsvector([
            playbook.label || '',
            playbook.description || '',
          ].join(' '))
          return { key: playbook.id, value: tsvector_intersect(search_tsvector, playbook_tsvector).size }
        }))
        playbooks.sort((a, b) => search_scores[b.id] - search_scores[a.id])
      }
      return playbooks
    }
  })
  .build()

export const PublicUserPlaybooks = API.get('/api/v1/public/user/playbooks')
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
    if ('db' in db && db.db instanceof PgDatabase) {
      const embedding = search ? await embed({ openai: await openai, content: search }) : undefined
      const playbooks = (await db.db.raw(subst => `
        select "published_playbook".*
        from "published_playbook"
        ${search ? `left join fpl_embedding on "published_playbook".playbook = fpl_embedding.id` : ''}
        ${inputs || outputs ? `where ${[
          inputs && `"published_playbook".inputs like any(${subst(inputs.map(s => `%${s}%`))})`,
          outputs && `"published_playbook".outputs like any(${subst(outputs.map(s => `%${s}%`))})`,
        ].filter(v => !!v).join(' and ')}` : ''}
        ${search ? `order by "fpl_embedding".embedding <-> ${subst(embedding)}` : `order by "clicks" desc`}
        offset ${subst(skip)}
        limit ${subst(limit)}
      `)).rows.map(row => schema.published_playbook.codec.decode(row))
      return playbooks
    } else {
      let playbooks = await db.objects.published_playbook.findMany({
        orderBy: {
          clicks: 'desc',
        },
        skip,
        take: limit,
      })
      if (inputs) playbooks = playbooks.filter(playbook => !inputs.some(spec => !playbook.inputs.split('\t').includes(spec)))
      if (outputs) playbooks = playbooks.filter(playbook => !outputs.some(spec => !playbook.outputs.split('\t').includes(spec)))
      if (search) {
        const search_tsvector = tsvector(search)
        const search_scores = dict.init(playbooks.map(playbook => {
          const playbook_tsvector = tsvector([
            playbook.label || '',
            playbook.description || '',
          ].join(' '))
          return { key: playbook.id, value: tsvector_intersect(search_tsvector, playbook_tsvector).size }
        }))
        playbooks.sort((a, b) => search_scores[b.id] - search_scores[a.id])
      }
      return playbooks
    }
  })
  .build()

export const UserPlaybooks = API.get('/api/v1/user/playbooks')
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
      skip: inputs.query.skip,
      take: inputs.query.limit,
    })
    return playbooks
  })
  .build()

export const UserPlaybook = API.get('/api/v1/user/playbooks/[id]')
  .query(z.object({
    id: z.string(),
  }))
  .call(async (inputs, req, res) => {
    const fpl = await fpprg.getFPL(inputs.query.id)
    if (fpl === undefined) throw new NotFoundError()
    // TODO: move this logic into a db trigger
    // TODO: address playbook republishing
    const [publicUserPlaybook] = await db.objects.user_playbook.findMany({
      where: {
        playbook: inputs.query.id,
        public: true,
      },
      take: 1,
      orderBy: {
        created: 'asc',
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
    if (!session || !session.user) return
    const userPlaybook = await db.objects.user_playbook.findUnique({
      where: {
        user: session.user.id,
        playbook: inputs.query.id,
      },
    })
    if (!userPlaybook) return
    return { public: userPlaybook.public }
  })
  .build()

export const UpdateUserPlaybook = API.post('/api/v1/user/playbooks/[id]/update')
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
    fpl = await fpprg.upsertFPL(fpl)
    if (fpl === undefined) throw new NotFoundError()
    fpl = await fpl.rebaseInTime()
    if (fpl === undefined) throw new NotFoundError()
    fpl = await fpprg.upsertFPL(fpl)
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const playbookInputs = fpl.resolve().map(cell => krg.getProcessNode(cell.process.type)).filter(node => 'prompt' in node && dict.isEmpty(node.inputs)).map(node => node.output.spec).join(', ')
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
    if (inputs.body.user_playbook.public) {
      // ensure published playbooks get an embedding computed
      await indexSoon({ db, fpl_id: fpl.id })
    }
    return fpl.id
  })
  .build()

export const PublishUserPlaybook = API.post('/api/v1/user/playbooks/[id]/publish')
  .query(z.object({
    id: z.string(),
  }))
  .body(z.object({
    public: z.boolean(),
  }))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    if (inputs.body.public) {
      await indexSoon({ db, fpl_id: inputs.query.id })
    }
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

export const DeleteUserPlaybook = API.post('/api/v1/user/playbooks/[id]/delete')
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
