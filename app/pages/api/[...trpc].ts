export { openApiHandler as default } from '@/app/openapi'

export const config = {
  api: {
    // workflows submitted to POST /api/db/fpl can embed (large) data
    bodyParser: { sizeLimit: '8mb' },
  },
}
