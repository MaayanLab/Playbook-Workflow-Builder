export function downloadUrl(url: string, filename?: string) {
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