import numpy as np
from voronoid import *
from drawvoronoid import *


# decorater used to block function printing to the console
def blockPrinting(func):
    def func_wrapper(*args, **kwargs):
        # block all printing to the console
        sys.stdout = open(os.devnull, 'w')
        # call the method in question
        value = func(*args, **kwargs)
        # enable all printing to the console
        sys.stdout = sys.__stdout__
        # pass the return value of the method back
        return value

    return func_wrapper


@blockPrinting
def findRegion(points, xmin=None, xmax=None, ymin=None, ymax=None):
    new_regions = []
    xmin, xmax, ymin, ymax, finalpoints, regions = voronoid(
        points, xmin, xmax, ymin, ymax)
    regions = np.asarray(regions)
    for region in regions:
        region = np.asarray(region)
        c = region.mean(axis=0)
        angles = np.arctan2(region[:, 1] - c[1], region[:, 0] - c[0])
        region = np.array(region)[np.argsort(angles)]
        #polygon = vertices[region]
        new_regions.append(region)
    print(new_regions)

    return new_regions, finalpoints, xmin, xmax, ymin, ymax


def plotVoronoi(points, xmin=None, xmax=None, ymin=None, ymax=None):
    regions, finalpoints, xmin, xmax, ymin, ymax = findRegion(
        points, xmin, xmax, ymin, ymax)
    for region in regions:
        plt.fill(*zip(*region), alpha=0.4)
        for ver in region:
            plt.annotate('({:.2f},{:.2f})'.format(
                ver[0], ver[1]), (ver[0], ver[1]))
    for p in finalpoints:
        plt.annotate('({:.2f},{:.2f})'.format(p[0], p[1]), (p[0], p[1]))

    plt.plot(finalpoints[:, 0], finalpoints[:, 1], 'ko')
    #for p in finalpoints:
    #    plt.annotate('({:.2f},{:.2f})'.format(p[0],p[1]), (p[0],p[1]))
    plt.xlim(xmin - 0.1, xmax + 0.1)
    plt.ylim(ymin - 0.1, ymax + 0.1)

    plt.show()

if __name__=='__main':
    np.random.seed(11)
    points = np.random.randint(0, 10, (5, 2))
    points = [(x[0], x[1]) for x in points]
    #points = [(0, 1), (1, 8), (9, 0), (9, 4), (10, 10)]
    #n = len(points)
    points = list(set(points))

    xygraph = Xygraph(vl=points)
    xmin = xygraph.xmin - 1
    ymin = xygraph.ymin - 1
    v = (xmin, ymin)
    points.sort(key=lambda p: ((p[0]-v[0])**2 +
                            (p[1]-v[1])**2)**(1/2), reverse=True)
    cur_points = []

    for i in range(n):
        print(cur_points)
        cur_points.append(points.pop())
        plotVoronoi(cur_points, -1, 13, -1, 11)
    
