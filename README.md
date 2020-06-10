# Voronoi-diagram
Implementation of voronoi diagram with incremental algorithm
## Description of files
* [line_intersection.py](/line_intersection.py): Determine whether 2 lines intersect
* [bisector.py](/bisector.py): Define how to create the bisector and the functions realated to it such as finding intersection
* [xygraph.py](/xygraph.py): Create the empty frame for the initial double connected edgelist
* [dcel.py](/dcel.py): 
  * Double connected edge list implemetation
  * Update new points
* [voronoid.py](/voronoid.py): Initalize the double connected edgelist and add new points
* [drawvoronoid.py](/drawvoronoid.py): Draw the voronoi diagram
* [plot.ipynb](/plot.ipynb): Import and visualize the datapoints
* [hospital.ipynb](/hospital.ipynb): Extract latitudes and longitudes of hospitals from the map then visualize them with Voronoi

## How to Import and Run the Algorithm
To try out this algorithm with your data points, run
```python
from voronoid import Xygraph
from drawvoronoid import plotVoronoi

points = [(0, 1), (1, 8), (9, 0), (9, 4), (10, 10)]
plotVoronoi(points, -1, 13, -1, 11)
```
![image](https://github.com/khuyentran1401/Voronoi-diagram/blob/master/Screenshot%20from%202020-06-08%2023-31-46.png)
