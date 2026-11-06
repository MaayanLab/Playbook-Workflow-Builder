import { createOpenApiNextHandler, generateOpenApiDocument } from 'trpc-to-openapi'
import { appRouter } from '@/app/trpc/router'
import { createNextContext } from '@/app/trpc/context'
import packageJson from '@/package.json'

export const openApiHandler = createOpenApiNextHandler({
  router: appRouter,
  createContext: createNextContext,
  onError: ({ error, path }) => {
    if (error.code === 'INTERNAL_SERVER_ERROR') console.error(`[openapi] ${path}:`, error)
  },
})

export function openApiDocument() {
  return generateOpenApiDocument(appRouter, {
    title: 'Playbook Workflow Builder API',
    description: 'Public APIs for creating, inspecting & exporting workflows and accessing GA4GH DRS objects.',
    version: packageJson.version,
    baseUrl: `${(process.env.PUBLIC_URL || '').replace(/\/$/, '')}/api`,
    docsUrl: `${(process.env.PUBLIC_URL || '').replace(/\/$/, '')}/api/openapi.json`,
    tags: ['GA4GH DRS', 'Workflows'],
  })
}
