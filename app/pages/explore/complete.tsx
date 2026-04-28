import krg from '@/app/krg'
import { useExRouter } from '@/app/fragments/ex-router'
import React from 'react'

export default function Explore() {
  const router = useExRouter()
  React.useEffect(() => {
    router.push(`/explore/Start/${krg.getDataNodes().map(node=>node.spec).join('/')}`)
  }, [])
  return null
}