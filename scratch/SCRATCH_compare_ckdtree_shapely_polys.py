import numpy as np
import numpy.random as r

from scipy.spatial import cKDTree
from shapely.geometry import Point, Polygon


#create 100 random points in 100x100 space
pts = zip(r.random(100)*100, r.random(100)*100)


#create function to query ckdtree
def tree(pts, rad):
    t = cKDTree(pts, leafsize = 100)
    q = t.query(pts, k = 2, distance_upper_bound = rad)
    return(q) 


#create function make a circle of radius rad around a point pt
def make_circ(pt, rad, angs):
    pts = zip(cos(angs)*rad+pt[0], sin(angs)*rad+pt[1])
    poly = Polygon(pts)
    return(poly)

#create a function to loop over a list of points and check whether other points within its radius
#NOTE: NEED TO VECTORIZE, AND SET TO PRINT LIST OF ONLY THOSE INDIVIDS' IDS WITHIN THE CIRCLE, EXCLUDING SELF
def check_circs(pts, rad, angs):
    pairs = {}
    for pt in range(pts.shape[0]):
        poly = make_circ(pts[pt,], rad, angs)
        within = [Point(pts[i,]).within(poly) for i in range(pts.shape[1])]
        pairs[pt] = within
    return(pairs)




#COMPARE PERFORMANCE!

