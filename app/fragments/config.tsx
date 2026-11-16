'use client'

import React from 'react'

/**
 * NOTE: `process.env.NEXT_PUBLIC_*` mustn't be read directly (outside of @/app/runtime-config):
 *  next inlines those at build time, while they're configured at runtime (so that one build can
 *  be deployed with different urls). Use `useRuntimeConfig()` (or `runtimeConfig()`) instead.
 */
const fallbackRuntimeConfig: Record<string, string> = {
  NEXT_PUBLIC_URL: '',
  NEXT_PUBLIC_LANDING_PAGE: '/graph/extend',
}
const RuntimeConfigContext = React.createContext(fallbackRuntimeConfig)

// the configuration last provided, for code outside of react (e.g. utils/download.ts)
let currentRuntimeConfig = fallbackRuntimeConfig

/**
 * Provides the runtime configuration (see @/app/runtime-config), passed down from the root layout
 */
export function RuntimeConfig({ value, children }: React.PropsWithChildren<{ value: Record<string, string> }>) {
  const runtimeConfig = React.useMemo(() => ({ ...fallbackRuntimeConfig, ...value }), [value])
  currentRuntimeConfig = runtimeConfig
  return (
    <RuntimeConfigContext.Provider value={runtimeConfig}>
      {children}
    </RuntimeConfigContext.Provider>
  )
}

export function useRuntimeConfig() {
  return React.useContext(RuntimeConfigContext)
}

/** The runtime configuration, outside of react components (prefer `useRuntimeConfig` within them) */
export function runtimeConfig() {
  return currentRuntimeConfig
}
