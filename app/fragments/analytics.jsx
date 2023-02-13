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

function Matomo({ url, siteId }) {
  const router = useRouter()
  React.useEffect(() => {
    var _paq = window._paq = window._paq || [];
    /* tracker methods like "setCustomDimension" should be called before "trackPageView" */
    _paq.push(['trackPageView']);
    _paq.push(['enableLinkTracking']);
    _paq.push(['setTrackerUrl', url+'matomo.php']);
    _paq.push(['setSiteId', siteId]);
  }, [])
  React.useEffect(() => {
    if (typeof window._paq === 'undefined') return
    window._paq.push(['setCustomUrl', router.asPath])
    window._paq.push(['trackPageView'])
  }, [router.asPath])
  return <Script
    src={`${url}/matomo.js`}
    strategy="lazyOnload"
  />
}

export default function Analytics() {
  return (
    <>
      {process.env.NEXT_PUBLIC_GA_MEASUREMENT_ID ?
        <GA id={process.env.NEXT_PUBLIC_GA_MEASUREMENT_ID} />
        : null}
      {process.env.NEXT_PUBLIC_MATOMO_URL && process.env.NEXT_PUBLIC_MATOMO_SITE_ID ?
        <Matomo url={process.env.NEXT_PUBLIC_MATOMO_URL} siteId={process.env.NEXT_PUBLIC_MATOMO_SITE_ID} />
        : null}
    </>
  )
}
