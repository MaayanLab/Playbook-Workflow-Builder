import Script from "next/script"

export default function GA() {
  return (
    <Script
      src={`https://www.googletagmanager.com/gtag/js?id=${process.NEXT_PUBLIC_GA_MEASUREMENT_ID}`}
      strategy="lazyOnload"
      onLoad={() => {
        window.dataLayer = window.dataLayer || [];
        function gtag(){window.dataLayer.push(arguments);}
        gtag('js', new Date());
        gtag('config', process.NEXT_PUBLIC_GA_MEASUREMENT_ID);
      }}
    />
  )
}
