import OpenAI from 'openai'
import cache from '@/utils/global_cache'

export default cache('openai', async () => {
  if (!process.env.OPENAI_API_KEY) {
    console.warn(`OPENAI_API_KEY not defined`)
    throw new Error('Missing OPENAPI_KEY')
  }
  return new OpenAI({ apiKey: process.env.OPENAI_API_KEY })
})
