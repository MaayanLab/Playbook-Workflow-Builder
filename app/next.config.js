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

module.exports = {
  experimental: {
    externalDir: true,
  },
  async redirects() {
    return [
      {
        source: '/',
        destination: '/testing',
        permanent: false,
      },
    ]
  }
}