import {BoxSelectTool, BoxSelectToolView} from "models/tools/gestures/box_select_tool";
import {ColumnDataSource} from "models/sources/column_data_source";
import {PanEvent} from "core/ui_events";
import * as p from "core/properties";
import { bk_tool_icon_box_select } from "styles/icons";
import { SelectionGeometry } from "core/bokeh_events";

export class EnhanceToolView extends BoxSelectToolView {

  // Overriding the PanEvent handling native to BoxSelectTool
  // This is done so we can use all of the features of BoxSelectTool
  // without the gray overlay obscuring the non-selected parts
  _pan_end(ev: PanEvent): void {
    const { sx, sy } = ev;
    const curpoint = [sx, sy];
    const [sxlim, sylim] = this._compute_limits(curpoint as [number, number]);
    this._do_select(sxlim, sylim, true, this._select_mode(ev));
    this.model.overlay.update({ left: null, right: null, top: null, bottom: null });
    this._base_point = null;
    this.plot_view.push_state('box_select', { selection: this.plot_view.get_selection() });
    // This is the piece of code that removes the gray overlay on the non-selected parts
    this._clear();
  }

  // This excerpt of code doesn't actually change anything it overrides.
  // It was/is going to be intended to change the event that is returned when
  // the selection is made.
  //
  // The purpose for making a new event would have been to add a field that
  // would distinguish between whether the selection is made in the distance
  // difference matrix display or the queued covariance submatrix display.
  // Currently, the event that is returned (SelectionGeometry) is
  // indistinguishable when returned from a selection in either display
  _emit_selection_event(geometry: any, final = true) {
    const { x_scale, y_scale } = this.plot_view.frame; 
    let geometry_data; 
    switch (geometry.type) { 
        case "point": { 
            const { sx, sy } = geometry; 
            const x = x_scale.invert(sx); 
            const y = y_scale.invert(sy); 
            geometry_data = Object.assign(Object.assign({}, geometry), { x, y }); 
            break; 
        } 
        case "span": { 
            const { sx, sy } = geometry; 
            const x = x_scale.invert(sx); 
            const y = y_scale.invert(sy); 
            geometry_data = Object.assign(Object.assign({}, geometry), { x, y }); 
            break; 
        } 
        case "rect": { 
            const { sx0, sx1, sy0, sy1 } = geometry; 
            const [x0, x1] = x_scale.r_invert(sx0, sx1); 
            const [y0, y1] = y_scale.r_invert(sy0, sy1); 
            geometry_data = Object.assign(Object.assign({}, geometry), { x0, y0, x1, y1 }); 
            break; 
        } 
        case "poly": { 
            const { sx, sy } = geometry; 
            const x = x_scale.v_invert(sx); 
            const y = y_scale.v_invert(sy); 
            geometry_data = Object.assign(Object.assign({}, geometry), { x, y }); 
            break; 
        } 
    } 
    // A new event object would need to be made and then used to replace
    // SelectionGeometry here
    this.plot_model.trigger_event(new SelectionGeometry(geometry_data, final));
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
  properties: EnhanceTool.Props;
  __view_type__: EnhanceToolView;

  constructor(attrs?: Partial<EnhanceTool.Attrs>) {
    super(attrs);
      this.tool_name = "Enhance Select";
      this.icon = bk_tool_icon_box_select;
      this.event_type = "pan";
      this.default_order = 30;
  }

  static init_EnhanceTool(): void {
    this.prototype.default_view = EnhanceToolView;

    this.define<EnhanceTool.Props>({
      source: [ p.Instance ],
    })
  }
}

