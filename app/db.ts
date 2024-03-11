import create_database from "@/utils/orm"
import * as schema from "@/db"
import cache from '@/utils/global_cache'

export default cache('db', () => create_database({
  connectionString: process.env.DATABASE_URL,
  schema: {
    user: schema.user,
    account: schema.account,
    session: schema.session,
    verification_token: schema.verification_token,
    user_integrations: schema.user_integrations,
    proxy_session: schema.proxy_session,
    suggestion: schema.suggestion,
    data: schema.data,
    process: schema.process,
    process_input: schema.process_input,
    resolved: schema.resolved,
    process_complete: schema.process_complete,
    fpl: schema.fpl,
    cell_metadata: schema.cell_metadata,
    playbook_metadata: schema.playbook_metadata,
    upload: schema.upload,
    user_upload: schema.user_upload,
    user_upload_complete: schema.user_upload_complete,
    user_playbook: schema.user_playbook,
    chat_message: schema.chat_message,
  },
}))
