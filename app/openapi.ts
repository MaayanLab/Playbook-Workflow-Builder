import { createOpenApiFetchHandler, generateOpenApiDocument } from 'trpc-to-openapi'
import { appRouter } from '@/app/trpc/router'
import { createFetchContext } from '@/app/trpc/context'
import packageJson from '@/package.json'

/** Serve the OpenAPI procedures, at their paths within `endpoint` (e.g. `/api/bco/...`, `/` for the site's root) */
export function openApiFetchHandler(req: Request, endpoint: `/${string}` = '/api') {
  return createOpenApiFetchHandler({
    router: appRouter,
    req,
    endpoint,
    createContext: () => createFetchContext(req),
    onError: ({ error, path }) => {
      if (error.code === 'INTERNAL_SERVER_ERROR') console.error(`[openapi] ${path}:`, error)
    },
  })
}

export function openApiDocument() {
  return generateOpenApiDocument(appRouter, {
    title: 'Playbook Workflow Builder API',
    description: 'Public APIs for creating, inspecting & exporting workflows, uploading files and accessing GA4GH DRS objects.',
    version: packageJson.version,
    baseUrl: `${(process.env.PUBLIC_URL || '').replace(/\/$/, '')}/api`,
    docsUrl: `${(process.env.PUBLIC_URL || '').replace(/\/$/, '')}/api/openapi.json`,
    tags: ['GA4GH DRS', 'Workflows', 'Uploads'],
    // API keys are created on the account page (docs/user/api.md)
    securitySchemes: {
      apiKey: { type: 'apiKey', in: 'header', name: 'x-api-key' },
      bearer: { type: 'http', scheme: 'bearer', description: 'An API key, as `Authorization: Bearer <key>`' },
    },
  })
}
