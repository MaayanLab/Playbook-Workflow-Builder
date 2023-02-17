import React from 'react'
import Icon from '@/app/components/icon'
import { share_icon } from '@/icons'
import usePublicUrl from '@/utils/next-public-url'

export default function ShareButton({ id }: { id: string | undefined }) {
  const publicUrl = usePublicUrl({ absolute: true })
  const [share, setShare] = React.useState(false)
  const onClick = React.useCallback(() => {
    const graphUrl = document.getElementById('graph-url') as HTMLInputElement
    graphUrl.select()
    graphUrl.setSelectionRange(0, 99999)
    navigator.clipboard.writeText(graphUrl.value)
  }, [id])
  React.useEffect(() => { if (share) { onClick() } }, [share])
  if (!share) {
    return (
      <button className={`bp4-button bp4-minimal`} onClick={() => {setShare(true)}}>
        <Icon icon={share_icon} color="black" />
      </button>
    )
  } else {
    return (
      <div className={`bp4-control-group inline-block${share ? '': ' hidden'}`}>
        <input id="graph-url" type="text" className="bp4-input" value={`${publicUrl}/report${id ? `/${id}` : ''}`} readOnly />
        <button className="bp4-button bp4-icon-link" onClick={onClick} />
        <button className="bp4-button bp4-icon-cross" onClick={() => {setShare(false)}} />
      </div>
    )
  }
}