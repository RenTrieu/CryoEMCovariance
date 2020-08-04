from bokeh.plotting import figure, save
from bokeh.models import ColorBar, LinearColorMapper, Plot, Range1d, LinearAxis, FixedTicker, FuncTickFormatter

color_mapper = LinearColorMapper(palette="Spectral6", low=0, high=60)
ticker = FixedTicker(ticks=[5,15,25,35,45,55])
formatter = FuncTickFormatter(code="""
    data = {5: '0-10', 15: '10-20', 25: '20-30', 35: '30-40', 45: '40-50', 55: '50+'}
    return data[tick]
""")

cbar = ColorBar(color_mapper=color_mapper, ticker=ticker, formatter=formatter,
                major_tick_out=0, major_tick_in=0, major_label_text_align='left',
                major_label_text_font_size='10pt', label_standoff=2)

p = Plot(x_range=Range1d(0,1), y_range=Range1d(0,1), width=500, height=500, toolbar_location=None)
p.add_layout(LinearAxis(), 'left')
p.add_layout(LinearAxis(), 'below')

p.add_layout(cbar)

save(p)
