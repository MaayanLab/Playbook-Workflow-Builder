import { z } from 'zod'
import { authedProcedure, router } from '@/app/trpc/init'

/** Typed API owned by the Elysium component. */
const elysiumRouter = router({
  signUpload: authedProcedure
    .input(z.object({ uid: z.string() }))
    .mutation(async ({ input }) => {
      if (!process.env.ELYSIUM_USERNAME || !process.env.ELYSIUM_PASSWORD) throw new Error('Not Configured')
      const search = new URLSearchParams({ username: process.env.ELYSIUM_USERNAME, password: process.env.ELYSIUM_PASSWORD })
      const response = await fetch(`https://maayanlab.cloud/charon/signpolicy?${search}`)
      if (!response.ok) throw new Error(await response.text())
      const signed = await response.json() as Record<string, string>
      return {
        url: `https://${signed.bucket}.s3.amazonaws.com/`,
        formData: {
          key: `${signed.uid}/${input.uid}-\${filename}`,
          AWSAccessKeyId: signed.cid,
          acl: 'private',
          success_action_redirect: 'success.html',
          policy: signed.policy,
          signature: signed.signature,
          'Content-Type': 'application/octet-stream',
        },
      }
    }),
})

export default elysiumRouter
