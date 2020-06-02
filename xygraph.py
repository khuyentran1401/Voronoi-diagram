"""
xygraph implements a 2D map formed by undirected edges between
vertices.
"""

__author__ = "Angel Yanguas-Gil"

import math as m


class Xygraph:
    """Represents a set of vertices connected by undirected edges.
    The vertices are stored in a list of coordinates, while
    the edges are stored as a pair of indices (i,j) of the vertices
    list.
    """

    def __init__(self, vl=[], el=[]):
        """Creates the 2D graph formed by a list of vertices (x,y)
        and a list of indices (i,j)
        """
        self.vl = vl
        self.el = el
        if self.vl != []:
            self.minmax()

    def minmax(self):
        """Determines the boundary box of the vertices in the graph"""
        vx = [v[0] for v in self.vl]
        vy = [v[1] for v in self.vl]
        self.xmax, self.xmin = max(vx), min(vx)
        self.ymax, self.ymin = max(vy), min(vy)
