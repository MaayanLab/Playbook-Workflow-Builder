import { openApiFetchHandler } from '@/app/openapi'

export const dynamic = 'force-dynamic'
export const maxDuration = 300

// the OpenAPI handler routes by the request url (it doesn't use the route params); wrapped, as
//  route handlers are called with (request, context), not the handler's (request, endpoint)
function handler(request: Request) {
  return openApiFetchHandler(request, '/api')
}

export { handler as GET, handler as POST, handler as PUT, handler as PATCH, handler as DELETE }
