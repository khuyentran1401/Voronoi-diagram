# Voronoi-diagram
Implementation of voronoi diagram with incremental algorithm
## Description of files
* [bisector.py](/bisector.py): Define how to create the bisector and the functions realated to it
* [xygraph.py](/xygraph.py): Create the empty frame for the initial dcel
* [dcel.py](/dcel.py): 
  * Double connected edge list implemetation
  * Update new points
* [voronoid.py](/voronoid.py): Initalize the double connected edgelist and add new points
* [drawvoronoid.py](/drawvoronoid.py): Draw the voronoi diagram
* [plot.ipynb](/plot.ipynb): Import and visualize the datapoints

## How to Import and Run the Algorithm
To try out this algorithm with your data points, run
```python
from voronoid import Xygraph
from drawvoronoid import plotVoronoi

points = [(0, 1), (1, 8), (9, 0), (9, 4), (10, 10)]
plotVoronoi(points, -1, 13, -1, 11)
```
