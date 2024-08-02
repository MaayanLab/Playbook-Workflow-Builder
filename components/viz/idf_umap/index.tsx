import python from '@/utils/python'
import { PlotlyPlot } from '@/components/viz/plotly'
import { MetaNode } from '@/spec/metanode'
import { plot_icon } from '@/icons'
import { GMT } from '@/components/data/gene_matrix_transpose'
import { DMT } from '@/components/data/drug_matrix_transpose'

export const UMAPFromXMT = [GMT, DMT].map(XMT =>
  MetaNode(`UMAPFrom[${XMT.spec}]`)
    .meta({
      label: `UMAP plot from ${XMT.meta.label}`,
      description: `Construct UMAP plot with ${XMT.meta.label}`,
      icon: [plot_icon],
    })
    .inputs({ matrix: XMT })
    .output(PlotlyPlot)
    .resolve(async (props) => await python(
      'components.viz.idf_umap.idf_umap',
      { kargs: [props.inputs.matrix] },
      message => props.notify({ type: 'info', message }),
    ))
    .story(props => ({
      abstract: `UMAP was applied to the inverse document frequency of the ${XMT.meta.label.toLocaleLowerCase()}.`,
      introduction: `The Inverse Document Fequency (IDF) is an adjustment for the fact that some genes appear more frequently then others, and is used commonly in TF-IDF for information retrieval of words\\ref{doi:10.1017/CBO9781139058452.002}.  Uniform Manifold Approximation and Projection (UMAP) is a dimensionality reducation technique widely used for visualization of high dimensional data\\ref{doi:10.48550/arXiv.1802.03426}.`,
      methods: `The ${XMT.meta.label.toLocaleLowerCase()} in ${props.input_refs?.matrix} is visualized by first computing IDF\\ref{doi:10.1017/CBO9781139058452.002} followed by UMAP\\ref{doi:10.48550/arXiv.1802.03426}.`,
      legend: `A scatter plot of the ${XMT.meta.label.toLocaleLowerCase()} in ${props.input_refs?.matrix} visualized by IDF\\ref{doi:10.1017/CBO9781139058452.002} followed by UMAP\\ref{doi:10.48550/arXiv.1802.03426}.`,
    }))
    .build()
)
