import React from 'react'
import * as Bokeh from '@bokeh/bokehjs'
import useAsyncEffect from 'use-async-effect'

/**
 * Use bokeh embed_item on a div with a unique id.
 * Note that embed_item appends the plot to the div, BokehPlot (below)
 *  thus replaces the component entirely when the plot prop changes.
 */
function BokehEmbeddedItem({ plot }: { plot: any }) {
    const id = React.useId()
    useAsyncEffect(async () => {
        if (!plot) return
        try {
            // Give the element to bokeh to manage
            const views = await Bokeh.embed.embed_item(plot, id)
            return () => {
                // on unmount, we tell bokeh to stop managing this element
                if (views) views.forEach(item => item.remove())
            }
        } catch (e) {
            console.warn(e)
        }
    }, [plot])
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
