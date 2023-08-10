export function downloadUrl(url: string, filename?: string) {
  if (url.startsWith((process.env.NEXT_PUBLIC_URL || '').replace(/https?:/, 'drs:'))) {
    url = url.replace(/^drs:\/\/([^\/]+)\/(.+)$/, `${process.env.NEXT_PUBLIC_URL||''}/ga4gh/drs/v1/objects/$2/access/https/data`)
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