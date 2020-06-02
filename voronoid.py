from bisector import *
from line_intersection import *
from dcel import *
from xygraph import *
from sklearn.neighbors import NearestNeighbors
import numpy as np
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.patches as patches
from PIL import Image
import numpy as np
import matplotlib 
from PIL import Image, ImageDraw

#matplotlib.use('GTK3Cairo')



def voronoid(points):

    xygraph = Xygraph(vl=points)
    xmin = xygraph.xmin - 1
    xmax = xygraph.xmax + 3
    ymin = xygraph.ymin - 1
    ymax = xygraph.ymax + 1


    # Ensure there is no points in the same location
    points = list(set(points))
    n = len(points)

    # Sort points in advance by the distance to the l
    v = (xmin, ymin)

    points.sort(key=lambda p: ((p[0]-v[0])**2 +
                               (p[1]-v[1])**2)**(1/2), reverse=True)

    # Points displayed
    cur_points = []

    p = points.pop()

    print(p)
    cur_points.append(p)
    cur_points = cur_points
    print('border', xmin, xmax, ymin, ymax)

    V = Dcel(vl=[(xmin, ymax), (xmin, ymin), (xmax, ymin), (xmax, ymax)], el=[
             (0, 1), (1, 2), (2, 3), (3, 0)], site=p, border=[xmin, xmax, ymin, ymax])

   

    # change this
    for i in range(n-1):

        print('\n')
        print('i', i)

        p = points.pop()

        print(p)
        cur_points.append(p)
        # Closest site to the new site p_i+1
        nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(
            np.array(cur_points))
        indices = nbrs.kneighbors(
            np.array(p).reshape(1, -1), return_distance=False)
        pc = cur_points[indices[0][1]]

        # face associated with pc
        fn = V.getFace(pc)

        p1, q1 = perpendicular_bisector(p, pc, xmin, xmax, ymin, ymax)
        print('point closest', pc)

        intersect_vl = []
        intersect_edges = {}

        vertices = [v.coord for v in V.vertices]

        move = -10e-2
        m = {}
        findEdge = True
        eps = 10e-3

        for h in fn.hedges:

            m[h] = slope(h.v1.coord, h.origin.coord)

            if doIntersect(p1, q1, h.vertices[0].coord, h.vertices[1].coord):

                # Find the intersection between the bisector and the intersect line

                pt = intersection(
                    p1, q1, h.vertices[0].coord, h.vertices[1].coord)
                
                # Handle the intersect vertex and is the same as the existing vertex
                if pt in vertices and findEdge:
                    shift = 0.0000000001
                    pt_before = pt
                    pt = intersection(
                        p1, q1, h.vertices[0].coord, h.vertices[1].coord, shift)
                
                    if abs(pt_before[0]-pt[0]) >  eps:
                        if pt_before[0]-pt[0] > 0:
                            pt = (-eps, pt[1])
                        else:
                            pt = (eps, pt[1])
                    
                   
                    if not isOnLine(pt, h):
                        continue
                    else:
                        findEdge = False


                vertex = Vertex(pt[0], pt[1])

                intersect_vl.append(vertex)
                print('intersect vertex', vertex)
                intersect_edges[vertex] = h
                print('intersect hedge', h)

        V.update(p, pc, intersect_vl, intersect_edges, xmin, xmax, ymin, ymax)
    vertices = np.array([v.coord for v in V.vertices])
    #faces = np.array([f.site for f in V.faces])
    reg = [f.hedges for f in V.faces.values()]
    hedgelist = list(set([h for heds in reg for h in heds]))
    hedges = [(V.vertices.index(h.vertices[0]),
               V.vertices.index(h.vertices[1])) for h in hedgelist]
    

    #Take one pair of hedge
    ridge_vertices = []
    for h in hedges:
        if (h[1], h[0]) not in ridge_vertices:
            ridge_vertices.append(h)

    ridge_points = []
    for h in hedgelist:
        site = []
        if h.newface:
            site.append(cur_points.index(h.newface.site))
        else:
            site.append(-1)
        if h.twin.newface:
            site.append(cur_points.index(h.twin.newface.site))
        else:
            site.append(-1)
        
        ridge_points.append(site)
        
    regions = []
    for f in V.faces.values():
        region = []
        for v in f.vertices:            
            try:
                region.append(V.vertices.index(v))
            except:
                print(v.coord)
                pass
        
        regions.append(list(set(region)))

    index_region = [site for site in V.faces.keys()]
    point_region = [index_region.index(p) for p in cur_points]

    

    print(vertices)

    print(ridge_vertices)

    print(regions)

    print(ridge_points)

    print(point_region)

    return xmin, xmax, ymin, ymax, np.array(cur_points), vertices, np.array(ridge_points), np.array(ridge_vertices), np.array(regions), np.array(point_region)

    
    
    
    


if __name__ == '__main__':
    points = [(-5,-4),(-6, 4), (4, 3), (6,2)]
    # Handle the vertex intersection
    #points = [(-3, 1), (-1, 3), (1,3)]
    #np.random.seed(0)
#
    #x = np.random.randint(-10, 10, 3)
    #y = np.random.randint(-10, 10, 3)

    
    #points = list(zip(x,y))
    print(points)

    # Handle point on voronoi line
    #points= [(0, 2), (4, 2), (2, 2)]
    voronoid(points)

