import type { PropsWithChildren } from 'react'
import Shell from '@/app/app/_clients/shell'

/**
 * The site's shell (header, scrolling page, chat sidebar & the playbook state) for all the pages
 *  within: being a layout, it persists while navigating between them (rather than each page
 *  rendering its own, which remounted it all, and lost its state, on every navigation).
 */
export default function ShellLayout({ children }: PropsWithChildren) {
  return <Shell>{children}</Shell>
}
