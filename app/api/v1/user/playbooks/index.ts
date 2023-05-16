import { z } from 'zod'
import { API } from '@/spec/api'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { NotFoundError, UnauthorizedError } from '@/spec/error'
import db from '@/app/db'
import fpprg from '@/app/fpprg'

export const CreateUserPlaybook = API('/api/v1/user/playbooks/[id]')
  .query(z.object({
    id: z.string(),
    public: z.boolean().optional().transform(v => v === true),
  }))
  .body(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const fpl = await fpprg.getFPL(inputs.query.id)
    if (fpl === undefined) throw new NotFoundError()
    const ret = await db.objects.user_playbook.upsert({
      where: {
        user: session.user.id,
        playbook: inputs.query.id,
      },
      create: {
        user: session.user.id,
        playbook: inputs.query.id,
        title: fpl.metadata?.label,
        description: fpl.metadata?.description,
        public: inputs.query.public,
      },
      update: {
        title: fpl.metadata?.label,
        description: fpl.metadata?.description,
        public: inputs.query.public,
      },
    })
    if (!ret) throw new NotFoundError()
    return ret
  })
  .build()

export const UpdateUserPlaybook = API('/api/v1/user/playbooks/[id]/update/[new_id]')
  .query(z.object({
    id: z.string(),
    new_id: z.string(),
    public: z.boolean().optional().transform(v => v === true),
  }))
  .body(z.any())
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    
    const currentUserPlaybook = await db.objects.user_playbook.findUnique({
      where: {
        user: session.user.id,
        playbook: inputs.query.id,
      },
    })
    if (!currentUserPlaybook) throw new NotFoundError()
    const fpl = await fpprg.getFPL(inputs.query.new_id)
    if (fpl === undefined) throw new NotFoundError()
    const ret = await db.objects.user_playbook.update({
      where: {
        user: session.user.id,
        playbook: inputs.query.new_id,
      },
      data: {
        playbook: inputs.query.new_id,
        title: fpl.metadata?.label,
        description: fpl.metadata?.description,
        public: inputs.query.public,
      },
    })
    if (!ret) throw new NotFoundError()
    return ret
  })
  .build()

export const DeleteUserPlaybook = API('/api/v1/user/playbooks/[id]/delete')
  .query(z.object({
    id: z.string(),
  }))
  .body(z.any())
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user) throw new UnauthorizedError()
    const ret = await db.objects.user_playbook.delete({
      where: {
        user: session.user.id,
        playbook: inputs.query.id,
      },
    })
    if (!ret) throw new NotFoundError()
    return ret
  })
  .build()
