import React from 'react'
import useSWR from 'swr/immutable'
import fetcher from '@/utils/next-rest-fetcher'
import { useRouter } from 'next/router'

const fallbackRuntimeConfig = {
  NEXT_PUBLIC_URL: process.env.NEXT_PUBLIC_URL ? process.env.NEXT_PUBLIC_URL : '',
  NEXT_PUBLIC_LANDING_PAGE: process.env.NEXT_PUBLIC_LANDING_PAGE ? process.env.NEXT_PUBLIC_LANDING_PAGE : '/graph/extend',
}
const RuntimeConfigContext = React.createContext(fallbackRuntimeConfig)

export function RuntimeConfig({ children }: React.PropsWithChildren<{}>) {
  const router = useRouter()
  const { data: runtimeConfig } = useSWR<typeof fallbackRuntimeConfig>(() => router.isReady ? `${router.basePath}/api/config` : undefined, fetcher)
  return (
    <RuntimeConfigContext.Provider value={runtimeConfig ? runtimeConfig : fallbackRuntimeConfig}>
      {children}
    </RuntimeConfigContext.Provider>
  )
}

export function useRuntimeConfig() {
  return React.useContext(RuntimeConfigContext)
}
