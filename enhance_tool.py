from bokeh.core.properties import Instance
from bokeh.io import output_file, show
from bokeh.models import ColumnDataSource, Tool
from bokeh.plotting import figure
from bokeh.util.compiler import TypeScript
import random

output_file('tool.html')

TS_CODE = """
import {BoxSelectTool, BoxSelectToolView} from "models/tools/gestures/box_select_tool"
import {ColumnDataSource} from "models/sources/column_data_source"
import {PanEvent} from "core/ui_events"
import * as p from "core/properties"

export class EnhanceToolView extends BoxSelectToolView {

  // Overriding the PanEvent handling native to BoxSelectTool
  _pan_end(ev: PanEvent): void {
    const { sx, sy } = ev;
    const curpoint = [sx, sy];
    const [sxlim, sylim] = this._compute_limits(curpoint as [number, number]);
    this._do_select(sxlim, sylim, true, this._select_mode(ev));
    this.model.overlay.update({ left: null, right: null, top: null, bottom: null });
    this._base_point = null;
    this.plot_view.push_state('box_select', { selection: this.plot_view.get_selection() });
    this._clear();
  }
}

export namespace EnhanceTool {
  export type Attrs = p.AttrsOf<Props>

  export type Props = BoxSelectTool.Props & {
    source: p.Property<ColumnDataSource>
  }
}

export interface EnhanceTool extends EnhanceTool.Attrs {}

export class EnhanceTool extends BoxSelectTool {
  properties: EnhanceTool.Props
  __view_type__: EnhanceToolView

  constructor(attrs?: Partial<EnhanceTool.Attrs>) {
    super(attrs)
  }

  tool_name = "Enhance Span"
  icon = "bk-tool-icon-box_select"
  event_type = "pan" as "pan"
  default_order = 30

  static init_EnhanceTool(): void {
    this.prototype.default_view = EnhanceToolView

    this.define<EnhanceTool.Props>({
      source: [ p.Instance ],
    })
  }
}

"""


class EnhanceTool(Tool):
    __implementation__ = TypeScript(TS_CODE)
    source = Instance(ColumnDataSource)


source = ColumnDataSource(data=dict(x=[], y=[]))

plot = figure(x_range=(0, 10), y_range=(0, 10), tools=[EnhanceTool(source=source)])
plot.title.text = "Drag to draw on the plot"
#plot.line('x', 'y', source=source)
x = [random.randint(0, 10) for i in range(0,10)]
y = [random.randint(0, 10) for i in range(0,10)]
plot.scatter(x, y)

show(plot)
