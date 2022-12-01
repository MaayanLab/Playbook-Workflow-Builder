from bokeh.plotting import figure, show
from bokeh.transform import factor_cmap
from bokeh.models import ColumnDataSource, Select, NumeralTickFormatter
from bokeh.models.callbacks import CustomJS
from bokeh.layouts import column, row
import numpy as np
import seaborn as sns

def generate_colors(input_df, feature):
    pal = sns.color_palette()
    color = factor_cmap(feature, palette=pal.as_hex(), factors=input_df[feature].unique().tolist())

    return color

def interactive_circle_plot(input_df, x_lab, y_lab, feature):
    input_df['legend'] = input_df[feature]
    source = ColumnDataSource(input_df)

    TOOLTIPS = [
        ("index", "$index"),
        ("(x,y)", "($x, $y)"),
        ("feature", '@'+feature),
    ]

    n = input_df[feature].nunique()
    p = figure(height=400,tooltips=TOOLTIPS,x_axis_label=x_lab, y_axis_label=y_lab,sizing_mode="scale_width")
    point_size = 10 if input_df.shape[0] < 100 else 5
    color = generate_colors(input_df, feature)
    p1 = p.circle('x', 'y', size=point_size, source=source, legend_field='legend',fill_color= color, line_color=color)
    p.add_layout(p.legend[0], 'right')

    p.output_backend = "svg"
    p.xgrid.visible = False
    p.ygrid.visible = False
    p.xaxis.minor_tick_line_color = None
    p.yaxis.minor_tick_line_color = None
    p.xaxis.axis_label_text_font_size = '12pt'
    p.yaxis.axis_label_text_font_size = '12pt'
    p.xaxis[0].formatter = NumeralTickFormatter(format="0.0")
    p.yaxis[0].formatter = NumeralTickFormatter(format="0.0")

    return generate_plot(p, p1, source, input_df)

def generate_plot(p, p_circle, source_data, input_df):
    column_options = list(input_df.columns[:-4])

    my_colors = {}
    for feature in column_options:
        my_colors[feature] = generate_colors(input_df, feature)

    arg = dict(colors=my_colors, plot=p_circle, source = source_data,
           data_category = input_df.reset_index().to_dict(orient='list'), feature = input_df.columns[-2])

    select = Select(title="Plot to show:", options=column_options, width=800)
    select.js_on_change("value", CustomJS(args=arg, code="""

    plot.glyph.line_color = colors[this.value]
    plot.glyph.fill_color = colors[this.value]
    source.data['legend'] = data_category[this.value]
    source.data[feature] = data_category[this.value]

    """))

    layout = column(select, p)
    layout.sizing_mode = "stretch_both"
    return layout

