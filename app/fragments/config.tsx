import React from 'react'
import { useQuery } from '@tanstack/react-query'
import { useTRPC } from '@/app/trpc/client'

const fallbackRuntimeConfig = {
  NEXT_PUBLIC_URL: process.env.NEXT_PUBLIC_URL ? process.env.NEXT_PUBLIC_URL : '',
  NEXT_PUBLIC_LANDING_PAGE: process.env.NEXT_PUBLIC_LANDING_PAGE ? process.env.NEXT_PUBLIC_LANDING_PAGE : '/graph/extend',
}
const RuntimeConfigContext = React.createContext(fallbackRuntimeConfig)

export function RuntimeConfig({ children }: React.PropsWithChildren<{}>) {
  const trpc = useTRPC()
  const { data: runtimeConfig } = useQuery({ ...trpc.system.runtimeConfig.queryOptions(), staleTime: Infinity })
  return (
    <RuntimeConfigContext.Provider value={{ ...fallbackRuntimeConfig, ...runtimeConfig }}>
      {children}
    </RuntimeConfigContext.Provider>
  )
}

export function useRuntimeConfig() {
  return React.useContext(RuntimeConfigContext)
}
