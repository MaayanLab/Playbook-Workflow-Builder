import React from 'react'
import { useQuery } from '@tanstack/react-query'
import krg from '@/app/krg'
import { useTRPCClient } from '@/app/trpc/client'
import { withSuggestions } from '@/app/extensions/suggestions'

export const suggestionsQueryKey = (session_id?: string) => ['catalog-suggestions', session_id] as const

/**
 * The KRG: the components along with all user suggestions
 *  (so that playbooks containing suggested steps can be displayed)
 */
export default function useKRG({ session_id }: { session_id?: string }) {
  const trpc = useTRPCClient()
  const { data: suggestions } = useQuery({
    // A tunneled session has its own database and therefore its own suggestions.
    queryKey: suggestionsQueryKey(session_id),
    queryFn: () => trpc.catalog.suggestions.query(undefined, { context: { session_id } }),
    retry: false,
  })
  return React.useMemo(() => withSuggestions(krg, suggestions ?? []), [suggestions])
}
