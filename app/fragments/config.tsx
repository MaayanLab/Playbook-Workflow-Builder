'use client'

import React from 'react'

const fallbackRuntimeConfig = {
  NEXT_PUBLIC_URL: process.env.NEXT_PUBLIC_URL ? process.env.NEXT_PUBLIC_URL : '',
  NEXT_PUBLIC_LANDING_PAGE: process.env.NEXT_PUBLIC_LANDING_PAGE ? process.env.NEXT_PUBLIC_LANDING_PAGE : '/graph/extend',
}
const RuntimeConfigContext = React.createContext(fallbackRuntimeConfig)

/**
 * Provides the runtime configuration (see @/app/runtime-config), passed down from the root layout
 */
export function RuntimeConfig({ value, children }: React.PropsWithChildren<{ value: Record<string, string> }>) {
  const runtimeConfig = React.useMemo(() => ({ ...fallbackRuntimeConfig, ...value }), [value])
  return (
    <RuntimeConfigContext.Provider value={runtimeConfig}>
      {children}
    </RuntimeConfigContext.Provider>
  )
}

export function useRuntimeConfig() {
  return React.useContext(RuntimeConfigContext)
}
