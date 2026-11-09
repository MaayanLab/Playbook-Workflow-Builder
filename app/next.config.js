const fs = require('fs')
const path = require('path')
const dotenv = require('dotenv')
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
const env = dotenv.parse(fs.readFileSync(path.join(root, '.env')))
env.PYTHON_ROOT = env.PYTHON_ROOT || root
for (const key in env) {
  if (!(key in process.env)) {
    process.env[key] = env[key]
  }
}

// sane setup for Better Auth's public URL
if (!process.env.PUBLIC_URL) process.env.PUBLIC_URL = 'http://127.0.0.1:3000'
if (!process.env.NEXT_PUBLIC_URL) process.env.NEXT_PUBLIC_URL = process.env.PUBLIC_URL
if (!process.env.BETTER_AUTH_URL) process.env.BETTER_AUTH_URL = process.env.PUBLIC_URL
if (!process.env.LANDING_PAGE) process.env.LANDING_PAGE = '/graph/extend'
if (!process.env.NEXT_PUBLIC_LANDING_PAGE) process.env.NEXT_PUBLIC_LANDING_PAGE = process.env.LANDING_PAGE

/** @type {import('next').NextConfig} */
module.exports = {
  // a separate build directory lets a local CAVATICA worker (`next dev`) run alongside the app
  distDir: process.env.NEXT_DIST_DIR || '.next',
  images: {
    unoptimized: true,
  },
  experimental: {
    externalDir: true,
  },
  typescript: {
    ignoreBuildErrors: process.env.STRICT_TS!=='true',
  },
  async redirects() {
    return [
      {
        source: '/',
        destination: process.env.LANDING_PAGE ?? '/graph/extend',
        permanent: false,
      },
    ]
  },
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
  webpack: (config, { isServer, webpack }) => {
    if (isServer) {
      // database drivers (incl. native modules) are loaded from node_modules at runtime
      //  rather than bundled, files outside of app/ (externalDir) are otherwise bundled
      config.externals = [
        ...[].concat(config.externals ?? []),
        ...['pg', 'pg-boss', 'pg-listen', 'better-sqlite3', 'kysely', 'kysely/migration']
          .map(pkg => ({ [pkg]: `commonjs ${pkg}` })),
      ]
    } else {
      // server-only modules which end up in the client module graph through
      //  getServerSideProps re-exports, they're never used client-side
      config.resolve.alias = {
        ...config.resolve.alias,
        'better-sqlite3': false,
        'pg-boss': false,
      }
      config.resolve.fallback = {
        ...config.resolve.fallback,
        'pg-native': false,
        'util/types': false,
        tls: false,
        dns: false,
        net: false,
        fs: false,
        child_process: false,
        stream: false,
      };
    }
    return config
  }
}
