'use client'

import { useParams, usePathname, useRouter as useAppRouter, useSearchParams } from 'next/navigation'
import Link from 'next/link'
import React from "react";

/**
 * A router with the (pages router style) API the fragments were written against,
 *  implemented with the app router: `query` merges the route params & search params,
 *  `push` with `shallow` keeps the current search params.
 */
export function useExRouter() {
  const appRouter = useAppRouter()
  const pathname = usePathname()
  const params = useParams()
  const searchParams = useSearchParams()
  const query = React.useMemo(() => {
    const query: Record<string, string | string[]> = {}
    for (const [key, value] of (searchParams ?? new URLSearchParams()).entries()) query[key] = value
    for (const [key, value] of Object.entries(params ?? {})) {
      if (typeof value === 'string' || Array.isArray(value)) query[key] = value
    }
    return query
  }, [params, searchParams])
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
    appRouter.push(newUrl.pathname + newUrl.search + newUrl.hash, { scroll: options?.scroll })
    return Promise.resolve(true)
  }, [appRouter])
  const replace = React.useCallback((url: string) => {
    appRouter.replace(url)
    return Promise.resolve(true)
  }, [appRouter])
  return {
    pathname: pathname ?? '',
    asPath: `${pathname ?? ''}${searchParams?.size ? `?${searchParams}` : ''}`,
    query,
    isReady: true,
    push,
    replace,
    back: appRouter.back,
  }
}

export type ExRouter = ReturnType<typeof useExRouter>

/**
 * A link: use it (styled as a button where it's meant to look like one, rather than nesting a
 *  <button> in it) so that it behaves as one, e.g. shift/ctrl/middle click opens it elsewhere.
 * `disabled` links don't navigate (anchors have no disabled attribute).
 */
export function ExLink({ href, shallow, disabled, className, onClick, ...props }: React.ComponentProps<typeof Link> & { href: string, disabled?: boolean }) {
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
  if (disabled) {
    return <Link
      href={newUrl}
      {...props}
      aria-disabled
      tabIndex={-1}
      className={`${className ?? ''} bp6-disabled btn-disabled pointer-events-none`.trim()}
      onClick={evt => evt.preventDefault()}
    />
  }
  return <Link href={newUrl} className={className} onClick={onClick} {...props} />
}
