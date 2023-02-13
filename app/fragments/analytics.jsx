import React from 'react'
import { useRouter } from "next/router";
import Script from "next/script"

function GA({ id }) {
  React.useEffect(() => {
    window.dataLayer = window.dataLayer || [];
    function gtag(){window.dataLayer.push(arguments);}
    gtag('js', new Date());
    gtag('config', id);
  }, [])
  React.useEffect(() => {
    if (typeof window.dataLayer === 'undefined') return
    function gtag(){window.dataLayer.push(arguments);}
    gtag({
      event: 'pageview',
      pageUrl: router.asPath,
    })
  }, [router.asPath])
  return <Script
    src={`https://www.googletagmanager.com/gtag/js?id=${id}`}
    strategy="lazyOnload"
  />
}

export default function Analytics() {
  return (
    <>
      {process.env.NEXT_PUBLIC_GA_MEASUREMENT_ID ?
        <GA id={process.env.NEXT_PUBLIC_GA_MEASUREMENT_ID} />
        : null}
    </>
  )
}
