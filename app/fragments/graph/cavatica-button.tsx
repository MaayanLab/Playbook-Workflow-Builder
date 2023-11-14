import React from 'react'
import dynamic from 'next/dynamic'
import { connection_icon } from '@/icons'
import { useSessionWithId } from '@/app/extensions/next-auth/hooks'
import { useAPIMutation, useAPIQuery } from '@/core/api/client'
import { UserIntegrationsCAVATICA, UserIntegrationsCAVATICADisconnect, UserIntegrationsCAVATICALaunch } from '@/app/api/client'
import { useRouter } from 'next/router'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function CAVATICAButton({ session_id }: { session_id?: string }) {
  const router = useRouter()
  const userSession = useSessionWithId()
  const { data: userIntegrations } = useAPIQuery(UserIntegrationsCAVATICA, {}, { shouldRetryOnError: false })
  const { trigger: launchCAVATICA } = useAPIMutation(UserIntegrationsCAVATICALaunch, {})
  const { trigger: disconnectCAVATICA } = useAPIMutation(UserIntegrationsCAVATICADisconnect, { session_id })
  return (
    <>
      <button
        className="bp5-button bp5-minimal z-20"
        disabled={!userSession.data?.user}
        onClick={async () => {
          if (!userIntegrations?.cavatica_api_key) {
            router.push('/account/cavatica', '/account/cavatica', { shallow: true })
          } else if (!session_id) {
            const new_session_id = await launchCAVATICA()
            router.push(`/session/[session_id]/graph/extend`, `/session/${new_session_id}/graph/extend`, { shallow: true })
          } else {
            await disconnectCAVATICA()
            router.push(`/graph/extend`, `/graph/extend`, { shallow: true })
          }
        }}
      >
        <Icon
          icon={connection_icon}
          title={
            !userSession.data?.user ? 'Sign in to Connect with CAVATICA'
            : session_id ? 'Close Session'
            : !userIntegrations?.cavatica_api_key ? 'Configure CAVATICA Integration'
            : 'Connect to CAVATICA'
          }
          className={
            !userSession.data?.user ? 'fill-gray-400'
            : !userIntegrations?.cavatica_api_key ? 'fill-black dark:fill-white'
            : !session_id ? 'fill-black dark:fill-white'
            : 'fill-red-600'
          }
        />
      </button>
    </>
  )
}
