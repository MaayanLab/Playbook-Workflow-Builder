import React from 'react'
import usePublicUrl from '@/utils/next-public-url'
//  NOTE: the version of @bokeh/bokehjs must match the python bokeh (requirements.txt)
import bokehjsPackage from '@bokeh/bokehjs/package.json'

const bokehVersion: string = bokehjsPackage.version

/**
 * BokehJS is loaded from its own prebuilt bundles (served by @/app/app/bokehjs) rather than
 *  bundled with the app: BokehJS identifies its models by their class names, which the app's
 *  production build minifies away (e.g. "could not resolve type 'DocumentConfig'").
 */
type BokehEmbedded = { roots: { remove(): void }[] }
type BokehGlobal = { embed: { embed_item(item: unknown, target_id: string): Promise<BokehEmbedded> } }
declare global { interface Window { Bokeh?: BokehGlobal } }

const loaded = new Map<string, Promise<void>>()
function loadScript(src: string) {
  let promise = loaded.get(src)
  if (promise === undefined) {
    promise = new Promise<void>((resolve, reject) => {
      const script = document.createElement('script')
      script.src = src
      script.async = false
      script.onload = () => resolve()
      script.onerror = () => { loaded.delete(src); reject(new Error(`Failed to load ${src}`)) }
      document.head.appendChild(script)
    })
    loaded.set(src, promise)
  }
  return promise
}

/** The bundles a plot needs: the main one first, the others extend it */
function bundlesFor(plot: unknown) {
  const json = JSON.stringify(plot)
  return [
    'bokeh.min.js',
    'bokeh-widgets.min.js',
    ...json.includes('"TableColumn"') ? ['bokeh-tables.min.js'] : [],
    ...json.includes('"webgl"') ? ['bokeh-gl.min.js'] : [],
  ]
}

async function loadBokeh(publicUrl: string, plot: unknown) {
  for (const bundle of bundlesFor(plot)) {
    await loadScript(`${publicUrl}/bokehjs/${bokehVersion}/${bundle}`)
  }
  if (!window.Bokeh) throw new Error('BokehJS did not load')
  return window.Bokeh
}

/**
 * Use bokeh embed_item on a div with a unique id.
 * Note that embed_item appends the plot to the div, BokehPlot (below)
 *  thus replaces the component entirely when the plot prop changes.
 */
function BokehEmbeddedItem({ plot }: { plot: any }) {
    const id = React.useId()
    const publicUrl = usePublicUrl()
    const [error, setError] = React.useState<string | undefined>()
    React.useEffect(() => {
        if (!plot) return
        let cancelled = false
        let views: BokehEmbedded | undefined
        const remove = () => { if (views) for (const view of views.roots) view.remove() }
        // Give the element to bokeh to manage
        loadBokeh(publicUrl, plot)
            .then(Bokeh => cancelled ? undefined : Bokeh.embed.embed_item(plot, id))
            .then(embedded => {
                views = embedded
                // unmounted before the plot finished embedding
                if (cancelled) remove()
            })
            .catch(e => {
                console.warn(e)
                if (!cancelled) setError(String(e?.message ?? e))
            })
        return () => {
            cancelled = true
            // on unmount, we tell bokeh to stop managing this element
            remove()
        }
    }, [plot, id, publicUrl])
    const plotVersion = plot?.version ?? plot?.doc?.version
    return <>
        {error ? <div className="alert alert-error prose max-w-none">
            The plot could not be displayed ({error}).
            {plotVersion && plotVersion !== bokehVersion ? ` It was created with Bokeh ${plotVersion} while this page uses Bokeh ${bokehVersion}: recompute it to update it.` : null}
        </div> : null}
        <div id={id} />
    </>
}

/**
 * Render a BokehPlot json item
 */
export default function BokehPlot({ plot }: { plot: any }) {
    // a new key for each plot remounts BokehEmbeddedItem when the plot changes (only then), as
    //  Bokeh's embed_item appends to its element rather than replacing what's there
    const counter = React.useRef(0)
    const key = React.useMemo(() => ++counter.current, [plot])
    return plot ? <BokehEmbeddedItem key={key} plot={plot} /> : null
}
