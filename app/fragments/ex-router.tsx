import { useRouter } from "next/router";
import Link from 'next/link'
import React from "react";

export function useExRouter() {
  const router = useRouter()
  const push = React.useCallback((url: string, as_url?: string, options?: { shallow?: boolean, scroll?: boolean }) => {
    const qs = new URLSearchParams(window.location.search)
    const newUrl = new URL(as_url ?? url, window.location.origin)
    qs.entries()
      .filter(([key, _value]) => !newUrl.searchParams.has(key))
      .forEach(([key, value]) => {
        newUrl.searchParams.set(key, value)
      })
    if (as_url) {
      return router.push(url, newUrl.toString(), options)
    } else {
      return router.push(newUrl.toString(), undefined, options)
    }
  }, [router])
  return { ...router, push }
}

export function ExLink({ href, ...props }: React.ComponentProps<typeof Link> & { href: string }) {
  const router = useExRouter()
  const newUrl = React.useMemo(() => {
    const qs = new URLSearchParams(window.location.search)
    const newUrl = new URL(href, window.location.origin)
    qs.entries()
      .filter(([key, _value]) => !newUrl.searchParams.has(key))
      .forEach(([key, value]) => {
        newUrl.searchParams.set(key, value)
      })
    return newUrl.toString()
  }, [href, router.query])
  return <Link href={newUrl} {...props} />
}
