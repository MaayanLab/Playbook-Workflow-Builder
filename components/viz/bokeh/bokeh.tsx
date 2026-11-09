import React from 'react'
import * as Bokeh from '@bokeh/bokehjs'

/**
 * Use bokeh embed_item on a div with a unique id.
 * Note that embed_item appends the plot to the div, BokehPlot (below)
 *  thus replaces the component entirely when the plot prop changes.
 */
function BokehEmbeddedItem({ plot }: { plot: any }) {
    const id = React.useId()
    React.useEffect(() => {
        if (!plot) return
        let cancelled = false
        let views: Awaited<ReturnType<typeof Bokeh.embed.embed_item>> | undefined
        // Give the element to bokeh to manage
        Bokeh.embed.embed_item(plot, id)
            .then(embedded => {
                views = embedded
                // unmounted before the plot finished embedding
                if (cancelled && views) views.forEach(item => item.remove())
            })
            .catch(e => console.warn(e))
        return () => {
            cancelled = true
            // on unmount, we tell bokeh to stop managing this element
            if (views) views.forEach(item => item.remove())
        }
    }, [plot, id])
    return <div id={id} />
}

/**
 * Render a BokehPlot json item
 */
export default function BokehPlot({ plot }: { plot: any }) {
    // store the currently rendered plot in state
    const [visablePlot, setVisablePlot] = React.useState({ key: 0, plot })
    // if plot is updated, we'll update the key as well
    React.useEffect(() => setVisablePlot(({ key }) => ({ key: key + 1, plot })), [plot])
    // the updated key invalidates BokehEmbeddedItem forcing a remount, this avoids
    //  the problem caused by the fact that Bokeh's embed_item appends
    return visablePlot.plot ? <BokehEmbeddedItem key={visablePlot.key} plot={visablePlot.plot} /> : null
}
