#!/usr/bin/env python

"""
ugid_wx.py:

A small wxPython utility app to visualize pyugrids, etc.
"""

import wx


## import the installed version
from wx.lib.floatcanvas import NavCanvas, FloatCanvas


class DrawFrame(wx.Frame):

    """
    A frame used for the ugrid viewer

    """

    #some parameters for drawing:

    background_color = (200,200,200) # grey

    label_size = 16
    label_color = 'black'
    label_background_color = background_color

    node_color = 'black'

    face_color = 'cyan'
    face_edge_color = 'black'

    edge_color = 'red'
    
    def __init__(self, *args, **kwargs):
        wx.Frame.__init__(self, *args, **kwargs)

        self.CreateStatusBar()

        # Add the Canvas
        Canvas = NavCanvas.NavCanvas(self,-1,
                                     size = (500,500),
                                     ProjectionFun = None,
                                     Debug = 0,
                                     BackgroundColor = self.background_color,
                                     ).Canvas
        
        self.Canvas = Canvas

        FloatCanvas.EVT_MOTION(self.Canvas, self.OnMove ) 

        self.Show()
        Canvas.ZoomToBB()

    def Draw_UGRID(self, grid):
        """
        Draws a UGRID Object
        """
        
        self.Canvas.ClearAll()
        # add the elements:
        nodes = grid.nodes
        # add the elements:
        for i, f in enumerate(grid.faces):
            face = nodes[f]
            self.Canvas.AddPolygon(face, FillColor=self.face_color, LineColor=self.face_edge_color, LineWidth=2)
            mid = face.mean(axis=0)
            self.Canvas.AddText(`i`, mid, Size=self.label_size, Position='cc')
        
        # add the edges:
        for i, e in enumerate(grid.edges):
            edge = nodes[e]
            self.Canvas.AddLine(edge, LineColor=self.edge_color, LineWidth=3)
            mid = edge.mean(axis=0)
            self.Canvas.AddText(`i`,
                                mid,
                                Size=self.label_size,
                                Position='cc',
                                Color=self.label_color,
                                BackgroundColor=self.label_background_color)
            
        # add the Nodes
        for i, n in enumerate(nodes):
            self.Canvas.AddText(`i`, n, Size=self.label_size, BackgroundColor=self.label_background_color)
        self.Canvas.AddPointSet(nodes, Diameter=5, Color=self.node_color) 
        self.Canvas.ZoomToBB()

    def OnMove(self, event):
        """
        Updates the status bar with the world coordinates

        """
        self.SetStatusText("%.2f, %.2f"%tuple(event.Coords))


if __name__ == "__main__":
    from pyugrid import test_examples

    app = wx.App(False)
    F = DrawFrame(None, title="UGRID Test App", size=(700,700) )
    #F.Draw_UGRID( test_examples.two_triangles() )
    F.Draw_UGRID( test_examples.twenty_one_triangles() )

    app.MainLoop()
    
    
    







