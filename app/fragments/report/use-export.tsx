'use client'
/**
 * Exports which can take a long time to produce (LaTeX bundles, report PDFs, ...) are queued
 *  on the server (`exports.create`) and their progress is followed over SSE (`exports.watch`),
 *  a long request would otherwise be dropped by the ingress. A toast reports their progress
 *  and offers the result once it's ready. Exports require signing in, the file is kept with
 *  the user's uploads.
 */
import React from 'react'
import { OverlayToaster, Spinner, type Toaster, type ToastProps } from '@blueprintjs/core'
import { useTRPCClient } from '@/app/trpc/client'
import usePublicUrl from '@/utils/next-public-url'
import type { ExportFormat, ExportState } from '@/app/exports'

let toaster: Promise<Toaster> | undefined
function getToaster() {
  if (toaster === undefined) toaster = OverlayToaster.create({ position: 'bottom-right', maxToasts: 5 })
  return toaster
}

export type ExportRequest = {
  fpl_id: string,
  format: ExportFormat,
  metadata?: { title?: string, description?: string },
  /** what's being exported, e.g. "Report PDF" */
  label: string,
  /** what to do once it's ready (default download) */
  then?: 'download' | 'overleaf',
}

function triggerDownload(href: string, filename: string) {
  const a = document.createElement('a')
  a.href = href
  a.download = filename
  a.rel = 'noopener'
  document.body.appendChild(a)
  a.click()
  a.remove()
}

export function useExport(session_id?: string) {
  const trpc = useTRPCClient()
  const publicUrl = usePublicUrl({ absolute: true })
  return React.useCallback(async ({ fpl_id, format, metadata, label, then = 'download' }: ExportRequest) => {
    const context = { session_id }
    const prefix = session_id ? `/api/socket/${session_id}` : ''
    const toaster = await getToaster()
    let subscription: { unsubscribe: () => void } | undefined
    const onDismiss = () => { subscription?.unsubscribe() }
    const key = toaster.show({ message: `Requesting ${label}…`, icon: <Spinner size={16} />, timeout: 0, onDismiss })
    const update = (props: Omit<ToastProps, 'onDismiss'>) => toaster.show({ ...props, onDismiss }, key)
    const fail = (error: string) => update({
      message: `${label} could not be produced: ${error}`,
      intent: 'danger',
      icon: 'error',
      timeout: 0,
    })
    try {
      const { id } = await trpc.exports.create.mutate({ fpl_id, format, metadata }, { context })
      subscription = trpc.exports.watch.subscribe({ id }, {
        context,
        onData: (state: ExportState) => {
          if (state.status === 'queued') {
            update({ message: `${label} is queued…`, icon: <Spinner size={16} />, timeout: 0 })
          } else if (state.status === 'running') {
            update({ message: `Preparing ${label}, this can take a few minutes…`, icon: <Spinner size={16} />, timeout: 0 })
          } else if (state.status === 'error') {
            subscription?.unsubscribe()
            fail(state.error ?? 'unknown error')
          } else if (state.status === 'complete' && !state.download) {
            // the user deleted the file (from their uploads) since
            subscription?.unsubscribe()
            fail('the file is no longer available')
          } else if (state.status === 'complete' && state.download) {
            subscription?.unsubscribe()
            const href = `${prefix}/api/ga4gh/drs/v1/objects/${state.download.id}/access/https/data`
            if (then === 'overleaf') {
              // overleaf fetches the bundle itself, so it needs an absolute url
              const snip_uri = `${publicUrl ?? ''}${href}`
              // opening a tab outside of a click would be blocked, so the user opens it
              update({
                message: `${label} is ready`,
                intent: 'success',
                icon: 'tick',
                timeout: 0,
                action: { text: 'Open in Overleaf', href: `https://www.overleaf.com/docs?snip_uri=${encodeURIComponent(snip_uri)}`, target: '_blank' },
              })
            } else {
              triggerDownload(href, state.download.filename)
              update({
                message: `${label} is ready, it's also kept with your uploads`,
                intent: 'success',
                icon: 'tick',
                timeout: 30_000,
                // in case the browser didn't start the download
                action: { text: 'Download', href, download: state.download.filename } as ToastProps['action'],
              })
            }
          }
        },
        onError: (error) => fail(error.message),
      })
    } catch (error) {
      fail(error instanceof Error ? error.message : String(error))
    }
  }, [trpc, session_id, publicUrl])
}
