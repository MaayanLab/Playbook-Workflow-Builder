import React from 'react'
import classNames from 'classnames'

export default function IFrame({ className, ...props }: React.DetailedHTMLProps<React.IframeHTMLAttributes<HTMLIFrameElement>, HTMLIFrameElement>) {

  const ref = React.useRef<HTMLIFrameElement>(null)
  const [loading, setLoading] = React.useState(true)
  React.useEffect(() => {
    if (!ref || !ref.current) return
    const onload = () => {
      setLoading(false)
    }
    ref.current.addEventListener('load', onload)
    return () => {
      if (ref.current) {
        ref.current.removeEventListener('load', onload)
      }
    }
  }, [ref])
  return (
    <>
      <progress className={classNames('progress w-full', { 'hidden': !loading })}></progress>
      <iframe
        ref={ref}
        className={classNames(className, { 'hidden': loading })}
        {...props}
      />
    </>
  )
}