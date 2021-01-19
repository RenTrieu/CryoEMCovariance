#!/usr/bin/env python3

from bokeh.events import SelectionGeometry

class SelectionGeometryPlus(SelectionGeometry):
    """ Essentially the same as SelectionGeometry, but also contains
        an extra String field 
    """
    event_name = 'selectiongeometryplus'

    def __init__(self, model, geometry=None, plus=None, final=True):
        self.geometry = geometry
        self.final = final
        self.plus = plus
        super().__init__(model=model)
