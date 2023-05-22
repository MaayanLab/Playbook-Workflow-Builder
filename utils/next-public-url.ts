import * as React from 'react'
import { useRouter } from "next/router"
import { useRuntimeConfig } from '@/app/fragments/config'

/**
 * Figure out the public url
 *
 * Absolute if you need the origin as well, otherwise it'll be the absolute path relative to the page root.
 */
export default function usePublicUrl({ absolute }: { absolute?: boolean } = {}) {
  const router = useRouter()
  const runtimeConfig = useRuntimeConfig()
  const { asPath: serverLocation, isReady } = router
  const publicUrl = React.useMemo(() => {
    if (absolute) {
      if (typeof window === 'undefined') {
        return runtimeConfig.NEXT_PUBLIC_URL
      } else {
        const clientLocation = isReady ? window.location.toString() : runtimeConfig.NEXT_PUBLIC_URL
        // scheme://origin/basePath[/app/path]?whatever=true => scheme://origin/basePath
        const [clientBasePath, ..._] = clientLocation.split(serverLocation)
        return clientBasePath
      }
    } else {
      return router.basePath
    }
  }, [absolute, serverLocation, isReady])
  return publicUrl
}
