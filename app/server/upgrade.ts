/**
 * Websocket upgrades are shared by several plugins on one http server:
 *  socket.io, tRPC subscriptions and next.js (HMR in development).
 *
 * Plugins claim the paths they handle; next.js only receives unclaimed upgrades
 *  (it would otherwise handle, and close, every upgrade request).
 */
import type { IncomingMessage } from 'http'

const claimed = new Set<string>()

/** Claim upgrade requests to `path` (and paths below it if it ends with `/`) */
export function claimUpgradePath(path: string) {
  claimed.add(path)
}

export function isClaimedUpgrade(req: IncomingMessage) {
  const pathname = new URL(req.url ?? '/', 'http://localhost').pathname
  for (const path of claimed) {
    if (path.endsWith('/') ? pathname.startsWith(path) : pathname === path) return true
  }
  return false
}
