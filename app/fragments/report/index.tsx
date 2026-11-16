import React from 'react'
import dynamic from 'next/dynamic'
import { usePathname } from 'next/navigation'
import { useKRG } from '@/app/fragments/metapath'
import { parsePathname } from '@/app/app/_lib/routes'

const Cells = dynamic(() => import('@/app/fragments/report/cells'))

/**
 * The report view of the route (rendered by the report layouts, so that it persists while
 *  navigating; input edits also only change the url, see report/prompt); the shell & playbook
 *  state are provided by @/app/app/(shell)/layout.
 */
export default function App() {
  const krg = useKRG()
  const pathname = usePathname() ?? ''
  const route = React.useMemo(() => parsePathname(pathname), [pathname])
  if (route.kind !== 'report') return null
  return (
    <main className="grow flex flex-col mx-4">
      {route.id ?
        <Cells
          session_id={route.session_id}
          krg={krg}
          id={route.id}
        />
        : <div className="alert alert-error">Page not found</div>}
    </main>
  )
}
