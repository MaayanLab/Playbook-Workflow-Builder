import krg from '@/app/krg'
import { useRouter } from 'next/router'
import React from 'react'

export default function Explore() {
  const router = useRouter()
  React.useEffect(() => {
    router.push(`/explore/Start/${krg.getDataNodes().map(node=>node.spec).join('/')}`)
  }, [])
  return null
}