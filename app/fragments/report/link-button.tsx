import React from 'react'
import { link_icon } from '@/icons'
import usePublicUrl from '@/utils/next-public-url'
import dynamic from 'next/dynamic'
import classNames from 'classnames'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function LinkButton({ id, disabled }: { id: string | undefined, disabled: boolean }) {
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
      <button className="bp4-button bp4-minimal" disabled={disabled} onClick={() => {setShare(true)}}>
        <Icon
          icon={link_icon}
          className={disabled ? 'fill-gray-400' : 'fill-black dark:fill-white'}
          title={disabled ? 'Save to Share Link' : 'Share Link'}
        />
      </button>
    )
  } else {
    return (
      <div className={classNames('bp4-control-group inline-block', { 'hidden': !share })}>
        <input id="graph-url" type="text" className="bp4-input" value={`${publicUrl}/report${id ? `/${id}` : ''}`} readOnly />
        <button className="bp4-button bp4-icon-clipboard" onClick={onClick} />
        <button className="bp4-button bp4-icon-cross" onClick={() => {setShare(false)}} />
      </div>
    )
  }
}