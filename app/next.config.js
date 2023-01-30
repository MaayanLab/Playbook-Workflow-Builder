const fs = require('fs')
const dotenv = require('dotenv')

// create .env from .env.example if not present
if (!fs.existsSync('../.env')) {
  fs.copyFileSync('../.env.example', '../.env')
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

module.exports = {
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
        destination: '/sitemap',
        permanent: false,
      },
    ]
  }
}