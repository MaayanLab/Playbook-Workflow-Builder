import React from 'react'
import { link_icon } from '@/icons'
import dynamic from 'next/dynamic'
import classNames from 'classnames'

const Icon = dynamic(() => import('@/app/components/icon'))

export default function ShareButton({ disabled = false }: { disabled?: boolean }) {
  const [share, setShare] = React.useState(false)
  const onClick = React.useCallback(() => {
    const graphUrl = document.getElementById('graph-url') as HTMLInputElement
    graphUrl.value = window.location.href
    graphUrl.select()
    graphUrl.setSelectionRange(0, 99999)
    navigator.clipboard.writeText(graphUrl.value)
  }, [])
  React.useEffect(() => { if (share) { onClick() } }, [share])
  if (!share) {
    return (
      <button className="bp5-button bp5-minimal" disabled={disabled} onClick={() => {setShare(true)}}>
        <Icon
          icon={link_icon}
          className={'fill-black dark:fill-white'}
          title={'Share Temporary Link'}
        />
      </button>
    )
  } else {
    return (
      <div className={classNames('bp5-control-group inline-block', { 'hidden': !share })}>
        <input id="graph-url" type="text" className="bp5-input" defaultValue="" readOnly />
        <button className="bp5-button bp5-icon-clipboard" onClick={onClick} />
        <button className="bp5-button bp5-icon-cross" onClick={() => {setShare(false)}} />
      </div>
    )
  }
}