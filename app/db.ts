import create_database from "@/utils/orm"
import * as schema from "@/db"
import cache from '@/utils/global_cache'

export default cache('db', create_database({
  connectionString: process.env.DATABASE_URL,
  schema: {
    user: schema.user,
    account: schema.account,
    session: schema.session,
    verification_token: schema.verification_token,
    suggestion: schema.suggestion,
  },
}))
