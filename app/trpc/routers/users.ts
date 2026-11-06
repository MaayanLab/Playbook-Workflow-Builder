import { TRPCError } from '@trpc/server'
import { z } from 'zod'
import db from '@/app/db'
import { selectUserUploadComplete } from '@/db/queries'
import { uuid } from '@/db'
import { authedProcedure, publicProcedure, router } from '../init'

const profileFields = ['image', 'name', 'email', 'affiliation'] as const

export const usersRouter = router({
  profile: publicProcedure
    .input(z.object({ id: z.string() }))
    .query(async ({ input }) => {
      const user = await db.selectFrom('user').select(profileFields).where('id', '=', input.id).executeTakeFirst()
      if (!user) throw new TRPCError({ code: 'NOT_FOUND' })
      return user
    }),

  me: authedProcedure
    .query(async ({ ctx }) => {
      const user = await db.selectFrom('user').select(profileFields).where('id', '=', ctx.session.user.id).executeTakeFirst()
      if (!user) throw new TRPCError({ code: 'NOT_FOUND' })
      return user
    }),

  update: authedProcedure
    .input(z.object({
      name: z.string().optional(),
      image: z.string().optional(),
      affiliation: z.string().optional(),
    }).strict())
    .mutation(async ({ ctx, input }) => {
      if (Object.keys(input).length === 0) return null
      return await db.updateTable('user')
        .set(input)
        .where('id', '=', ctx.session.user.id)
        .returning(profileFields)
        .executeTakeFirst() ?? null
    }),

  delete: authedProcedure
    .mutation(async ({ ctx }) => {
      const user = await db.deleteFrom('user').where('id', '=', ctx.session.user.id).returning('id').executeTakeFirst()
      if (!user) throw new TRPCError({ code: 'NOT_FOUND' })
      return null
    }),

  suggestions: authedProcedure
    .query(async ({ ctx }) => {
      return await db.selectFrom('suggestion').selectAll().where('user', '=', ctx.session.user.id).execute()
    }),

  createSuggestion: authedProcedure
    .input(z.object({
      name: z.string(),
      inputs: z.string(),
      output: z.string(),
      user: z.string(),
      description: z.string(),
    }))
    .mutation(async ({ ctx, input }) => {
      if (input.user !== ctx.session.user.id) throw new TRPCError({ code: 'UNAUTHORIZED' })
      await db.insertInto('suggestion').values({ ...input, id: uuid() }).execute()
      return null
    }),

  deleteSuggestion: authedProcedure
    .input(z.object({ id: z.string() }))
    .mutation(async ({ ctx, input }) => {
      const suggestion = await db.deleteFrom('suggestion')
        .where('id', '=', input.id)
        .where('user', '=', ctx.session.user.id)
        .returning('id')
        .executeTakeFirst()
      if (!suggestion) throw new TRPCError({ code: 'NOT_FOUND' })
      return null
    }),

  uploads: authedProcedure
    .query(async ({ ctx }) => {
      const uploads = await selectUserUploadComplete(db).where('user_upload.user', '=', ctx.session.user.id).execute()
      return uploads.map(upload => ({
        ...upload,
        url: `${process.env.NEXT_PUBLIC_URL?.replace(/https?:\/\//, 'drs://')}/${upload.id}`,
      }))
    }),

  deleteUpload: authedProcedure
    .input(z.object({ id: z.string() }))
    .mutation(async ({ ctx, input }) => {
      const upload = await db.deleteFrom('user_upload')
        .where('id', '=', input.id)
        .where('user', '=', ctx.session.user.id)
        .returning('id')
        .executeTakeFirst()
      if (!upload) throw new TRPCError({ code: 'NOT_FOUND' })
      return null
    }),
})
