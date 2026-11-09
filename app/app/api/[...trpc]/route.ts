import { openApiFetchHandler } from '@/app/openapi'

export const dynamic = 'force-dynamic'
export const maxDuration = 300

// the OpenAPI handler routes by the request url (it doesn't use the route params)
export { openApiFetchHandler as GET, openApiFetchHandler as POST, openApiFetchHandler as PUT, openApiFetchHandler as PATCH, openApiFetchHandler as DELETE }
