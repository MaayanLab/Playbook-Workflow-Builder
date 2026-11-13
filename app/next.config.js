const fs = require('fs')
const path = require('path')
const { PHASE_DEVELOPMENT_SERVER } = require('next/constants')
const { parseEnv } = require('util')
if (!process.env.APP_ROOT) process.env.APP_ROOT = path.dirname(__dirname)
const root = process.env.APP_ROOT

// create .env from .env.example if not present
if (!fs.existsSync(path.join(root, '.env'))) {
  const envExample = fs.readFileSync(path.join(root, '.env.example')).toString()
  // Auto-generate the Better Auth secret
  const crypto = require('crypto')
  fs.writeFileSync(path.join(root, '.env'), envExample.replace(
    /(\n)?#(BETTER_AUTH_SECRET=)(\r?\n)/,
    `$1$2${crypto.randomBytes(20).toString('hex')}$3`
  ))
}

// update environment with .env
const env = parseEnv(fs.readFileSync(path.join(root, '.env'), 'utf8'))
env.PYTHON_ROOT = env.PYTHON_ROOT || root
for (const key in env) {
  if (!(key in process.env)) {
    process.env[key] = env[key]
  }
}

// these are read at runtime (by the server, & passed to the browser by app/runtime-config.ts), so
//  that one build can be deployed at different urls: nothing may depend on them at build time
// sane setup for Better Auth's public URL
if (!process.env.PUBLIC_URL) process.env.PUBLIC_URL = 'http://localhost:3000'
if (!process.env.NEXT_PUBLIC_URL) process.env.NEXT_PUBLIC_URL = process.env.PUBLIC_URL
if (!process.env.BETTER_AUTH_URL) process.env.BETTER_AUTH_URL = process.env.PUBLIC_URL
// `/` itself redirects to the landing page, so it can't be the landing page
if (!process.env.LANDING_PAGE || process.env.LANDING_PAGE === '/') process.env.LANDING_PAGE = '/graph/extend'
if (!process.env.NEXT_PUBLIC_LANDING_PAGE || process.env.NEXT_PUBLIC_LANDING_PAGE === '/') process.env.NEXT_PUBLIC_LANDING_PAGE = process.env.LANDING_PAGE

/** @type {import('next').NextConfig} */
const config = {
  // a separate build directory lets a local CAVATICA worker (`next dev`) run alongside the app
  distDir: process.env.NEXT_DIST_DIR || '.next',
  images: {
    unoptimized: true,
  },
  typescript: {
    ignoreBuildErrors: process.env.STRICT_TS!=='true',
  },
  // (redirects are resolved when building, the landing page redirect is at runtime: app/app/page.tsx)
  async rewrites() {
    return [
      {
        source: '/ga4gh/:path*',
        destination: '/api/ga4gh/:path*', // The :path parameter isn't used here so will be automatically passed in the query
      },
      {
        source: '/mcp',
        destination: '/api/mcp', // The :path parameter isn't used here so will be automatically passed in the query
      },
    ]
  },
  // pg-listen's dependency pg-format reads a file at a path computed at runtime, which can't be bundled
  //  (pg & better-sqlite3 are external by default, see next/dist/lib/server-external-packages.jsonc)
  serverExternalPackages: ['pg-listen'],
  turbopack: {
    // the app imports code from outside of app/ (components/, core/, db/, ...)
    root,
    // Component packages intentionally co-locate their browser UI and server resolver
    //  in one MetaNode DSL definition: the browser build projects resolvers out of
    //  component entrypoints, the server build projects React views & prompts out.
    rules: Object.fromEntries(['*.ts', '*.tsx'].map(glob => [glob, [
      {
        condition: { all: ['browser', { path: /components[\\/].+\.tsx?$/ }] },
        loaders: [{ loader: path.join(root, 'cli', 'metanode-client-loader.cjs'), options: { target: 'client' } }],
      },
      {
        condition: { all: [{ not: 'browser' }, { path: /components[\\/].+\.tsx?$/ }] },
        loaders: [{ loader: path.join(root, 'cli', 'metanode-client-loader.cjs'), options: { target: 'server' } }],
      },
    ]])),
  },
}

module.exports = (phase) => phase === PHASE_DEVELOPMENT_SERVER ? {
  ...config,
  // `next dev` only trusts localhost for its dev resources (e.g. hot reloading): also trust the host
  //  the app is served at, e.g. an https tunnel (which the chatbot needs: OpenAI calls `${PUBLIC_URL}/mcp`)
  allowedDevOrigins: [new URL(process.env.PUBLIC_URL).hostname],
} : config
