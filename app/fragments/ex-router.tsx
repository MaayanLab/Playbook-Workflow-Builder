import { useRouter } from "next/router";
import Link from 'next/link'
import React from "react";

export function useExRouter() {
  const router = useRouter()
  const push = React.useCallback((url: string, as_url?: string, options?: { shallow?: boolean, scroll?: boolean }) => {
    const newUrl = new URL(as_url ?? url, window.location.origin)
    if (options?.shallow) {
      Array.from(new URLSearchParams(window.location.search).entries())
        .filter(([key, _value]) => !newUrl.searchParams.has(key))
        .forEach(([key, value]) => newUrl.searchParams.set(key, value))
      Array.from(newUrl.searchParams.entries())
        .filter(([_key, value]) => value === '')
        .forEach(([key, _value]) => newUrl.searchParams.delete(key))
    }
    if (as_url) {
      return router.push(url, newUrl.toString(), options)
    } else {
      return router.push(newUrl.toString(), undefined, options)
    }
  }, [router])
  return { ...router, push }
}

export function ExLink({ href, shallow, ...props }: React.ComponentProps<typeof Link> & { href: string }) {
  const router = useExRouter()
  const newUrl = React.useMemo(() => {
    const newUrl = new URL(href, window.location.origin)
    if (shallow) {
      Array.from(new URLSearchParams(window.location.search).entries())
        .filter(([key, _value]) => !newUrl.searchParams.has(key))
        .forEach(([key, value]) => newUrl.searchParams.set(key, value))
      Array.from(newUrl.searchParams.entries())
        .filter(([_key, value]) => value === '')
        .forEach(([key, _value]) => newUrl.searchParams.delete(key))
    }
    return newUrl.toString()
  }, [href, shallow, router.query])
  return <Link href={newUrl} {...props} />
}
