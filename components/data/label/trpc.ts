import { z } from 'zod'
import { publicProcedure, router } from '@/app/trpc/init'
import { metadataFromFile } from './server/metadata'
import { updateMetadataColumn } from './server/update-metadata'

const file = z.object({
  url: z.string(),
  size: z.number().optional(),
  sha256: z.string().optional(),
  filename: z.string(),
  // (null for files without one, as in the file codec, see components/core/file)
  description: z.string().nullish(),
})

/** Typed API owned by the label component. */
const labelRouter = router({
  metadata: publicProcedure
    .input(file)
    .query(async ({ input }) => await metadataFromFile(input)),

  update: publicProcedure
    .input(z.object({ file, data: z.record(z.string(), z.record(z.string(), z.string())) }))
    .mutation(async ({ input }) => await updateMetadataColumn(input)),
})

export default labelRouter
