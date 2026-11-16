import { handleMcp } from '@/app/mcp'
import { routeResponse } from '@/app/app/api/_lib/response'

function handler(request: Request) {
  return routeResponse(() => handleMcp(request))
}

export { handler as GET, handler as POST, handler as DELETE }
