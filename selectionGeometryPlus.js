import { SelectionGeometry } from "../../../core/bokeh_events";

let SelectionGeometryPlus = class SelectionGeometryPlus extends SelectionGeometry {
        constructor(geometry, matrixType, final) {
            super();
            this.geometry = geometry;
            this.matrixType = matrixType;
            this.final = final;
        }
        _to_json() {
            const { geometry, matrixType, final } = this;
            return Object.assign(Object.assign({}, super._to_json()), { geometry, final });
        }
};
SelectionGeometryPlus.__name__ = "SelectionGeometryPlus";
SelectionGeometryPlus = __decorate([
    event("selectiongeometryplus")
], SelectionGeometryPlus);
export { SelectionGeometryPlus };
