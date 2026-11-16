import type { PropsWithChildren } from 'react'
import Report from '@/app/app/_clients/report'
import ClientOnly from '@/app/fragments/client-only'

/**
 * The report view lives in this layout rather than in the page so that it persists while
 *  navigating (the page, a catch-all, is remounted when its path changes): the page only validates
 *  the route (see ./[[...path]]/page).
 */
export default function ReportLayout({ children }: PropsWithChildren) {
  return <><ClientOnly><Report /></ClientOnly>{children}</>
}
