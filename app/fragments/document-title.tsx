'use client'

import React from 'react'

/**
 * Set the document title from a client component (next/head is not supported by the app router)
 */
export default function DocumentTitle({ title }: { title: string }) {
  React.useEffect(() => {
    const previous = document.title
    document.title = title
    return () => { document.title = previous }
  }, [title])
  return null
}
