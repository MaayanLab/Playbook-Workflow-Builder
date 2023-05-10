const withBundleAnalyzer = require('@next/bundle-analyzer')({
  enabled: process.env.NEXT_ANALYZE === 'true',
})
const fs = require('fs')
const path = require('path')
const dotenv = require('dotenv')

// create .env from .env.example if not present
if (!fs.existsSync('../.env')) {
  const envExample = fs.readFileSync(path.join('..', '.env.example')).toString()
  // Auto-generate NEXTAUTH_SECRET
  const crypto = require('crypto')
  fs.writeFileSync(path.join('..', '.env'), envExample.replace(
    /(\n)?#(NEXTAUTH_SECRET=)(\r?\n)/,
    `$1$2${crypto.randomBytes(20).toString('hex')}$3`
  ))
}

// update environment with .env
const env = dotenv.parse(fs.readFileSync('../.env'))
env.PYTHON_ROOT = env.PYTHON_ROOT || '../'
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

module.exports = withBundleAnalyzer({
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
        destination: process.env.LANDING_PAGE,
        permanent: false,
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
})
