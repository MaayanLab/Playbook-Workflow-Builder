import { runtimeConfig } from '@/app/fragments/config'

export function downloadUrl(url: string, filename?: string) {
  const publicUrl = runtimeConfig().NEXT_PUBLIC_URL || ''
  if (url.startsWith('drs://localhost:3000') || url.startsWith(publicUrl.replace(/^https?:/, 'drs:'))) {
    const session_match = /^\/session\/([^\/]+)/.exec(window.location.pathname)
    url = url.replace(
      /^drs:\/\/([^\/]+)\/(.+)$/,
      `${publicUrl}${
        session_match !== null ? `/api/socket/${session_match[1]}`: ''
      }/ga4gh/drs/v1/objects/$2/access/https/data`
    )
  }
  const a = document.createElement('a')
  a.href = url
  a.download = filename ? filename : url.split('/').slice(-1)[0]
  a.setAttribute('style', 'visibility: hidden')
  document.body.appendChild(a)
  a.dispatchEvent(new MouseEvent('click', { bubbles: false, cancelable: false, view: window }))
  document.body.removeChild(a)
}

export function downloadBlob(blob: Blob, filename: string) {
  return downloadUrl(URL.createObjectURL(blob), filename)
}