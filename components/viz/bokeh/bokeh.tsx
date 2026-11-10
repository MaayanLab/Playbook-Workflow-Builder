import React from 'react'
// BokehJS's package entry point imports its modules by bare paths (which only resolve within
//  its own build), its modules import each other relatively: so the embedding API & the widgets
//  (registering their models, e.g. Select) are imported from their modules.
//  NOTE: the version of @bokeh/bokehjs must match the python bokeh (requirements.txt)
import { embed } from '@bokeh/bokehjs/build/js/lib/main'
import '@bokeh/bokehjs/build/js/lib/models/widgets/main'

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
        let views: Awaited<ReturnType<typeof embed.embed_item>> | undefined
        const remove = () => { if (views) for (const view of views.roots) view.remove() }
        // Give the element to bokeh to manage
        embed.embed_item(plot, id)
            .then(embedded => {
                views = embedded
                // unmounted before the plot finished embedding
                if (cancelled) remove()
            })
            .catch(e => console.warn(e))
        return () => {
            cancelled = true
            // on unmount, we tell bokeh to stop managing this element
            remove()
        }
    }, [plot, id])
    return <div id={id} />
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
