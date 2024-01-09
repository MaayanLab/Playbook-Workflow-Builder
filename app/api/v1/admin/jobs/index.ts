import { z } from 'zod'
import { API } from '@/spec/api'
import { getServerSessionWithId } from '@/app/extensions/next-auth/helpers'
import { UnauthorizedError } from '@/spec/error'
import db from '@/app/db'
import fpprg from '@/app/fpprg'
import { PgDatabase } from '@/utils/orm/pg'

export const AdminJobs = API.get('/api/v1/admin/jobs')
  .query(z.object({}))
  .call(async (inputs, req, res) => {
    const session = await getServerSessionWithId(req, res)
    if (!session || !session.user || session.user.id !== '00000000-0000-0000-0000-000000000000') throw new UnauthorizedError()
    const jobs = await (db as PgDatabase).jobs()
    return await Promise.all(jobs.map(async ({ data, ...job }) => {
      if (!data) return null
      const proc = await fpprg.getProcess(data.id)
      if (!proc) return null
      return { id: data.id, type: proc.type, ...job }
    }))
  })
  .build()
