const fs = require('fs')
const path = require('path')
const dotenv = require('dotenv')
if (!process.env.APP_ROOT) process.env.APP_ROOT = path.dirname(__dirname)
const root = process.env.APP_ROOT

// create .env from .env.example if not present
if (!fs.existsSync(path.join(root, '.env'))) {
  const envExample = fs.readFileSync(path.join(root, '.env.example')).toString()
  // Auto-generate NEXTAUTH_SECRET
  const crypto = require('crypto')
  fs.writeFileSync(path.join(root, '.env'), envExample.replace(
    /(\n)?#(NEXTAUTH_SECRET=)(\r?\n)/,
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

// sane setup for NEXTAUTH_URL to avoid redundancy
if (!process.env.PUBLIC_URL) process.env.PUBLIC_URL = 'http://127.0.0.1:3000'
if (!process.env.NEXT_PUBLIC_URL) process.env.NEXT_PUBLIC_URL = process.env.PUBLIC_URL
if (!process.env.NEXTAUTH_URL_INTERNAL) process.env.NEXTAUTH_URL_INTERNAL = 'http://127.0.0.1:3000'
if (!process.env.NEXTAUTH_URL) process.env.NEXTAUTH_URL = process.env.PUBLIC_URL
if (!process.env.LANDING_PAGE) process.env.LANDING_PAGE = '/graph/extend'
if (!process.env.NEXT_PUBLIC_LANDING_PAGE) process.env.NEXT_PUBLIC_LANDING_PAGE = process.env.LANDING_PAGE

/** @type {import('next').NextConfig} */
module.exports = {
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
    if (!isServer) {
      config.resolve.fallback = {
        ...config.resolve.fallback,
        'pg-native': false,
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
