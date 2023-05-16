import { z } from 'zod'
import type { CreateUserPlaybook as CreateUserPlaybook_ } from './v1/user/playbooks'
export const CreateUserPlaybook = { path: "/api/v1/user/playbooks/[id]", method: "POST" as const, call: z.custom<typeof CreateUserPlaybook_['call']>() }
import type { UpdateUserPlaybook as UpdateUserPlaybook_ } from './v1/user/playbooks'
export const UpdateUserPlaybook = { path: "/api/v1/user/playbooks/[id]/update/[new_id]", method: "POST" as const, call: z.custom<typeof UpdateUserPlaybook_['call']>() }
import type { DeleteUserPlaybook as DeleteUserPlaybook_ } from './v1/user/playbooks'
export const DeleteUserPlaybook = { path: "/api/v1/user/playbooks/[id]/delete", method: "POST" as const, call: z.custom<typeof DeleteUserPlaybook_['call']>() }