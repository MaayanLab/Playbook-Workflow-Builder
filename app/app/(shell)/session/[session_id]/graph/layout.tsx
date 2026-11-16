import type { PropsWithChildren } from 'react'
import Graph from '@/app/app/_clients/graph'
import ClientOnly from '@/app/fragments/client-only'

/**
 * The graph view lives in this layout rather than in the page so that it persists while navigating
 *  between graphs & nodes (the page, a catch-all, is remounted when its path changes): the page
 *  only validates the route (see ./[[...path]]/page).
 */
export default function GraphLayout({ children }: PropsWithChildren) {
  return <><ClientOnly><Graph /></ClientOnly>{children}</>
}
