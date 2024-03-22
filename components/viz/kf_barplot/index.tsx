import { PlotlyJson, PlotlyPlot } from '@/components/viz/plotly'
import { MetaNode } from '@/spec/metanode'
import { barchart_icon } from '@/icons'
import { TumorGeneExpression } from '@/components/service/kf'
import { ScoredTissues } from '@/components/core/input/scored'


export const BarplotFromKfExpression = 
  MetaNode(`BarplotFrom[${TumorGeneExpression.spec}]`)
    .meta({
      label: `Bar plot from ${TumorGeneExpression.meta.label}`,
      description: `Construct Bar plot with ${TumorGeneExpression.meta.label}`,
      icon: [barchart_icon],
    })
    .inputs({ terms: TumorGeneExpression, other_terms: ScoredTissues})
    .output(PlotlyPlot)
    .resolve(async (props) => {
        const tumor_data = await props.inputs.terms;
        const tissue_data = await props.inputs.other_terms;

        const tumor_names = tumor_data.map(({Disease}) => Disease);
        const zscores = tumor_data.map(({zscore})=>zscore);

        const tissue_terms = tissue_data.map(({term}) => term);
        const tissue_zscores = tissue_data.map(({zscore}) => zscore);
        

        const bar_plot_start_index = 0;
        const bar_plot_end_index = 10;

        const data : PlotlyJson['data'] = [
            {
                x: tumor_names.slice(bar_plot_start_index,bar_plot_end_index),
                y: zscores.slice(bar_plot_start_index,bar_plot_end_index),
               name: 'tumor',
               type: 'bar' 
            },
            {
                x: tissue_terms.slice(bar_plot_start_index,bar_plot_end_index),
                y: tissue_zscores.slice(bar_plot_start_index,bar_plot_end_index),
               name: 'tissues',
               type: 'bar' 
            }
        ]
        

        const layout: PlotlyJson['layout'] = {
            title: 'Gene Expression in Tumors and Tissues',
            barmode: 'group',
            xaxis: {title: 'Diseases'},
            yaxis: {title: 'zscores'},
            margin: { b: 200 },
        }

        return {data, layout};
    })
    .story(props =>
      `To visualize the level of expression across ${TumorGeneExpression.meta.label.toLocaleLowerCase()}, a bar plot was created${''/* [FIGURE]*/}.`
    )
    .build()