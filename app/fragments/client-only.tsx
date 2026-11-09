'use client'

import React from 'react'

/**
 * Render children in the browser only (they're written assuming a browser, e.g. `window`).
 *  Unlike `dynamic(..., { ssr: false })` this isn't a suspense boundary, so pages finish
 *  rendering on the server (e.g. `notFound()` sets a 404 status) before the response starts.
 */
export default function ClientOnly({ children }: React.PropsWithChildren<{}>) {
  const [mounted, setMounted] = React.useState(false)
  React.useEffect(() => { setMounted(true) }, [])
  return mounted ? <>{children}</> : null
}
