import { openApiFetchHandler } from '@/app/openapi'

export const dynamic = 'force-dynamic'

/**
 * GA4GH APIs (e.g. DRS, which resolves `drs://<host>/<id>` to `https://<host>/ga4gh/drs/v1/...`):
 *  the OpenAPI procedures whose paths start with /ga4gh, served from the site's root (they're also
 *  under /api, with the rest of the API). More specific routes (e.g. a DRS object's data, a
 *  download) take precedence over this one.
 */
function handler(request: Request) {
  return openApiFetchHandler(request, '/')
}

export { handler as GET, handler as POST }
