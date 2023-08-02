import React from 'react'
import { useAPIQuery } from '@/core/api/client'
import { UserIntegrationsCAVATICAStatus } from '../api/client'

export default function SessionStatus({ session_id }: { session_id: string }) {
  const { data: status, error } = useAPIQuery(UserIntegrationsCAVATICAStatus, { session_id }, { refreshInterval: 5 })
  if (!status) {
    if (!error) return <div className="alert alert-info prose prose-lg max-w-none justify-center">Waiting for status..</div>
    else return <div className="alert alert-error prose prose-lg max-w-none justify-center">{error.toString()}</div>
  }
  if (!status.state) return <div className="alert alert-info prose prose-lg max-w-none justify-center">Sending job to CAVATICA..</div>
  if (status.state === 'ERROR') return <div className="alert alert-error prose prose-lg max-w-none justify-center">An error occurred while submitting the task.</div>
  if (status.state === 'QUEUED') return <div className="alert alert-info prose prose-lg max-w-none justify-center">Task is queued..</div>
  if (status.state === 'INITIALIZING') return <div className="alert alert-info prose prose-lg max-w-none justify-center">CAVATICA is initializing the task.</div>
  if (status.state === 'RUNNING') return <div className="alert alert-info prose prose-lg max-w-none justify-center">CAVATICA has started the task, waiting to connect.</div>
  if (status.state === 'CONNECTED') return <div className="alert alert-success prose prose-lg max-w-none justify-center">All executions are now being sent to CAVATICA worker.</div>
  if (status.state === 'CANCELING') return <div className="alert alert-warning prose prose-lg max-w-none justify-center">Session is being cancelled.</div>
  if (status.state === 'CANCELED') return <div className="alert alert-error prose prose-lg max-w-none justify-center">Session was cancelled.</div>
  if (status.state === 'EXECUTOR_ERROR') return <div className="alert alert-error prose prose-lg max-w-none justify-center">An error occurred while executing.</div>
  if (status.state === 'COMPLETE') return <div className="alert alert-error prose prose-lg max-w-none justify-center">Session was closed.</div>
  return <div className="alert alert-error prose prose-lg max-w-none justify-center">Status is {status.state}.</div>
}
