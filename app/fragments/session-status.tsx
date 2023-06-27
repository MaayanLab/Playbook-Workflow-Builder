import React from 'react'

export default function SessionStatus({ session_id }: { session_id?: string }) {
  if (!session_id) return null
  return <div className="alert alert-success prose prose-lg max-w-none justify-center">All executions are being sent to CAVATICA worker.</div>
}
