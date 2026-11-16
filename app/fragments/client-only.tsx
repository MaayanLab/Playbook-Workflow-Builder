'use client'

import React from 'react'

const subscribe = () => () => {}

/**
 * Whether this renders in the browser: false on the server & while hydrating (so that hydration
 *  matches the server's html), true afterwards, and right away for what mounts after hydration
 *  (e.g. after navigating, so that there's no empty frame).
 */
export function useIsClient() {
  return React.useSyncExternalStore(subscribe, () => true, () => false)
}

/**
 * Render children in the browser only (they're written assuming a browser, e.g. `window`).
 *  Unlike `dynamic(..., { ssr: false })` this isn't a suspense boundary, so pages finish
 *  rendering on the server (e.g. `notFound()` sets a 404 status) before the response starts.
 */
export default function ClientOnly({ children }: React.PropsWithChildren<{}>) {
  return useIsClient() ? <>{children}</> : null
}
