# Copyright 2008, Angel Yanguas-Gil
from math import sqrt
from xygraph import Xygraph
from bisector import *
from line_intersection import *
import math as m
import numpy as np


class DcelError(Exception):
	pass


class Vertex:

	def __init__(self, x, y):
		self.x = x
		self.y = y
		self.coord = (x, y)
		self.hedgelist = []

	def sortincident(self):
		self.hedgelist.sort(key=lambda h: h.angle)

	def sortthree(self, new_site, close_site, h1=None):
		if h1:
			if h1 not in self.hedgelist:
				h1 = h1.twin
			print(h1)
		
		compare_hedge = self.hedgelist[:]
		if h1:
			compare_hedge.remove(h1)

		# Determines if a point is to the right of a hedge
		
		
		for i, h in enumerate(compare_hedge):
			site = siteClose(new_site, close_site, h)
			#if h.newface == None or h.newface.site == site:
				
			if  lefton(h, site):
				self.hedgelist.remove(h)
				self.hedgelist.insert(0, h)
			
				

	def __str__(self):
		return '({}, {})'.format(self.x, self.y)


class Hedge:
	"""Minimal implementation of a half-edge of a 2D dcel"""

	def __init__(self, v1, v2):
		# The .coord is defined as the vertex it points to
		self.v1 = v1
		self.origin = v2
		self.twin = None
		self.newface = None
		self.nexthedge = None
		self.angle = hangle(v2.x-v1.x, v2.y-v1.y)
		self.prevhedge = None
		self.length = m.sqrt((v2.x-v1.x)**2 + (v2.y-v1.y)**2)
		self.vertices = (v1, v2)

	def __str__(self):
		return 'Edge: {} -> {}'.format(self.v1, self.origin)


class Face:
	"""Implements a newface of a 2D dcel"""

	def __init__(self, site=None):
		self.wedge = None
		self.site = site
		self.external = None
		self.hedges = []
		self.vertices = []

	def __str__(self):
		return 'Face - site {}'.format(str(self.site))

	def area(self):
		h = self.wedge
		a = 0
		while(not h.nexthedge is self.wedge):
			p1 = h.origin
			p2 = h.nexthedge.origin
			a += p1.x*p2.y - p2.x*p1.y
			h = h.nexthedge

		p1 = h.origin
		p2 = self.wedge.origin
		a = (a + p1.x*p2.y - p2.x*p1.y)/2
		return a

	def vertexlist(self):
		h = self.wedge
		pl = [h.origin]
		while(not h.nexthedge is self.wedge):
			h = h.nexthedge
			pl.append(h.origin)
		return pl

	def isinside(self, p):
		"""Determines whether a point is inside a newface"""

		h = self.wedge
		inside = False
		if lefton(h, p):
			while(not h.nexthedge is self.wedge):
				h = h.nexthedge
				if not lefton(h, p):
					return False
			return True
		else:
			return False


class Dcel(Xygraph):
	"""
	Implements a doubly-connected edge list
	"""

	def __init__(self, border=[], vl=[], el=[], clip=None, site=None):
		Xygraph.__init__(self, vl, el)
		self.site = site
		self.bordervertices = []
		self.vertices = []
		self.hedges = []
		self.faces = {}
		if border != []:
			self.border = border
		if vl != []:
			self.build_dcel()

	def build_dcel(self):
		"""
		Creates the dcel from the list of vertices and edges
		"""
		

		# Step 1: vertex list creation
		for v in self.vl:
			self.vertices.append(Vertex(v[0], v[1]))

		# Step 2: hedge list creation. Assignment of twins and
		# vertices
		hedges_in = []
		for e in self.el:
			if e[0] >= 0 and e[1] >= 0:
				h1 = Hedge(self.vertices[e[0]], self.vertices[e[1]])
				h2 = Hedge(self.vertices[e[1]], self.vertices[e[0]])
				h1.twin = h2
				h2.twin = h1
				self.vertices[e[1]].hedgelist.append(h1)
				self.vertices[e[0]].hedgelist.append(h2)
				self.hedges.append(h2)
				self.hedges.append(h1)

				hedges_in.append(h1)

		# Step 3: Identification of next and prev hedges

		for v in self.vertices:
			# print('Vertex: ',v)
			# for edge in v.hedgelist:
				# print(edge)
			v.sortincident()

			# for edge in v.hedgelist:
			# print(edge)
			l = len(v.hedgelist)
			if l < 2:
				raise DcelError(
					"Badly formed dcel: less than two hedges in vertex")
			else:
				for i in range(l-1):
					v.hedgelist[i].nexthedge = v.hedgelist[i+1].twin
					v.hedgelist[i+1].twin.prevhedge = v.hedgelist[i]

				v.hedgelist[l-1].nexthedge = v.hedgelist[0].twin

				v.hedgelist[0].twin.prevhedge = v.hedgelist[l-1]

		# Step 4: Face assignment
		provlist = hedges_in[:]

		nf = 0
		nh = len(hedges_in)

		while nh > 0:
			h = provlist.pop()
			nh -= 1
			# We check if the hedge already points to a newface
			if h.newface == None:
				# print('h',h)
				f = Face(site=self.site)
				# print('h newface', f)

				nf += 1
				# We link the hedge to the new newface
				f.wedge = h
				f.hedges.append(h)
				f.wedge.newface = f
				f.vertices.append(h.origin)
				# And we traverse the boundary of the new newface
				while (not h.nexthedge is f.wedge):
					# print('cur:', h)
					h = h.nexthedge
					f.vertices.append(h.origin)
					f.hedges.append(h)
					# print('next:', h)
					h.newface = f
				self.faces[f.site] = f

		# And finally we have to determine the external newface
		for f in self.faces.values():
			f.external = f.area() < 0

	def update(self, new_site, close_site, intersect_vl, intersect_edges, xmin, xmax, ymin, ymax):


		ignore, update_hedges = self.two_points_update(new_site, close_site,
											   intersect_vl, intersect_edges)


		
		# If the intersect hedges are not the border, we will continue draw the bisector
		num_update = 1
		while update_hedges:
			# Check this
			num_update += 1
			for ver, hed in update_hedges.items():
				print(ver, hed)

				# Draw another bisector
				newface_close = hed.twin.newface
				print('\nClose face', newface_close)
				p = new_site
				pc = newface_close.site
				p1, q1 = perpendicular_bisector(p, pc, xmin, xmax, ymin, ymax)

				intersect_vl = []
				intersect_edges = {}
				findEdge = True
				eps = 10e-3

				for h in newface_close.hedges:
					if doIntersect(p1, q1, h.vertices[0].coord, h.vertices[1].coord):

						# Find the intersection between the bisector and the intersect line

						pt = intersection(
							p1, q1, h.vertices[0].coord, h.vertices[1].coord)
						if abs(pt[0] - ver.coord[0])/max(abs(pt[0]),1) < eps and abs(pt[1] - ver.coord[1])/max(abs(pt[1]),1) < eps:
							vertex = ver
							h = hed.twin

							intersect_vl.append(vertex)
							intersect_edges[vertex] = h

							# Have the first element of the list to be the intersected vertex
							if len(intersect_vl) == 2:
								second = intersect_vl.pop(0)
								intersect_vl.append(second)

						else:
							vertex = Vertex(pt[0], pt[1])
							intersect_vl.append(vertex)
							intersect_edges[vertex] = h

						print('intersection', vertex.coord)
						print('intersect hedge', h)

				ignore, update_hedges = self.two_points_update(p, pc,
													   intersect_vl, intersect_edges, ver, hed.twin, num_update, ignore)
				
	def two_points_update(self, new_site, close_site, intersect_vl, intersect_edges, intersected_ver=None, intersected_hed=None, num_update=None, ignoreface=False):
		update_vertices = []

		# Update vertex list
		for v in intersect_vl:	
			#Handle old vertex
			if v != intersected_ver:
				self.vertices.append(v)
				handle_vertex = v

			update_vertices.append(v)

		# Split the bisector into hb1 and hb2, append the hb1 and hb2 to 2 different vertices (intersection points)
		hb1 = Hedge(update_vertices[-2], update_vertices[-1])
		hb2 = Hedge(update_vertices[-1], update_vertices[-2])

		hb1.twin = hb2
		hb2.twin = hb1

		newface_hedges = [hb1, hb2]

		update_vertices[-1].hedgelist.append(hb1)
		update_vertices[-2].hedgelist.append(hb2)
		self.hedges += ([hb1, hb2])

		# For each intersect edge, split it into h1 and h2, link each for each vertex, delete the old edge

		new_hedges = []
		delete_hedges = []
		head, tail = None, None
		newface_hedges = []
		merge_hedge1 = None
		deletever = None
		mergevertex = False
		for i, v in enumerate(intersect_vl):
			print(i, v)

			if intersected_hed and intersect_edges[v] == intersected_hed:
	
				print('intersected hed', intersected_hed)
				to_compare = []
				for h in v.hedgelist:
					if h.v1 == intersected_hed.v1 or h.v1 == intersected_hed.origin:
						to_compare.append(h)
						print('to compare', h)
		
				deletehedge = todelete_hedge(
					to_compare[0], to_compare[1], handle_vertex.coord)
				delete_hedges.append(deletehedge)

				# Link between the last hedge to the new cut hedge
				for h in to_compare:
					if h != deletehedge:
						intersected_hed.prevhedge.nexthedge = h 
						h.prevhedge = intersected_hed.prevhedge
						print('cur', intersected_hed.prevhedge)
						print('next', intersected_hed.prevhedge.nexthedge)


				for ver in deletehedge.vertices:
					if ver != intersected_ver:
						print('delete ver', ver)
						deletever = ver 
						
					
					if ver != v:
						try:
							self.vertices.remove(ver)
						except:
							pass
						

				print('delete', deletehedge)
				print('test', deletehedge.twin.nexthedge)
				

				if ignoreface and i == 0:
					ignoreface = False
					mergevertex = True
					continue
				

				if isborder(self.border, intersect_edges[intersect_vl[1]]) or isborder(self.border, deletehedge.twin.nexthedge) or isborder(self.border, intersected_hed.nexthedge):
			
					if isborder(self.border, deletehedge.twin.nexthedge) or isborder(self.border, deletehedge.nexthedge):
						print('head', deletehedge.twin.nexthedge)
						print('tail', deletehedge.prevhedge)
						head = deletehedge.twin.nexthedge
						tail = deletehedge.prevhedge
						print('head1', intersected_hed.nexthedge.nexthedge)
					
					elif num_update >= 3 and isborder(self.border, intersect_edges[intersect_vl[1]]):
						print('iam here')
						change_edge = intersected_hed.nexthedge
						head = change_edge.nexthedge 
						tail = change_edge.twin.prevhedge
						print('head', head)
						print('tail',tail)

					merge_hedge1 = Hedge(tail.v1, head.origin)
					merge_hedge2 = Hedge(head.origin, tail.v1)
					merge_hedge1.twin = merge_hedge2
					merge_hedge2.twin = merge_hedge1
					print('new hedge', merge_hedge1)
					print('twin', merge_hedge2)

					head.origin.hedgelist.remove(head)
					head.origin.hedgelist.append(merge_hedge1)
					tail.v1.hedgelist.remove(tail.twin)
					tail.v1.hedgelist.append(merge_hedge1.twin)

					head.nexthedge.prevhedge = merge_hedge1
					merge_hedge1.nexthedge = head.nexthedge

					print('cur', merge_hedge1)
					print('next', merge_hedge1.nexthedge)

					tail.prevhedge.nexthedge = merge_hedge1
					merge_hedge1.prevhedge = tail.prevhedge
					print('prev', merge_hedge1.prevhedge)

					if num_update >= 3 and isborder(self.border, intersect_edges[intersect_vl[1]]):
						intersect_edges[intersect_vl[1]] = merge_hedge1
					
				else:
					ignoreface = True
				

			else:
				# Append the new merge point to the intersect if the old hedge is on the intersect hedge
				if merge_hedge1:
					if head == intersect_edges[v] or tail == intersect_edges[v]:
						intersect_edges[v] = merge_hedge1

				htail1, htail2, horigin1, horigin2 = split_hedge(
					v, intersect_edges[v])
				print('intersect', intersect_edges[v])

				# Update the twins of new splitting edges
				htail1.twin = htail2
				htail2.twin = htail1
				horigin1.twin = horigin2
				horigin2.twin = horigin1
				print('htail1', htail1, 'horigin1', horigin1)

				if head and tail and isOnLine(head.v1.coord, htail1):
					print('update link')
					htail1.nexthedge = head.nexthedge 
					head.nexthedge.prevhedge = htail1 
					print('cur', htail1)
					print('next', htail1.nexthedge)

					horigin1.prevhedge = tail.prevhedge
					tail.prevhedge.nexthedge = horigin1
					print('cur', tail.prevhedge)
					print('next', tail.prevhedge.nexthedge)

				if htail1.origin not in update_vertices and htail1.origin != deletever:
					update_vertices.append(htail1.origin)

				if horigin1.v1 not in update_vertices and horigin1.v1 != deletever:
					update_vertices.append(horigin1.v1)

				new_hedges += ([htail1, htail2, horigin1, horigin2])

				if not newface_hedges:
					newface_hedges = [htail1, horigin1]

		delete_hedges += ([h for h in intersect_edges.values()
						   if h != intersected_hed])
		delete_twins = [h.twin for h in delete_hedges]
		delete_hedges += delete_twins

		for i, v in enumerate(update_vertices):

			# Remove the split hedges in each old vertex
			if i >= 2 or v == intersected_ver:
				for h in delete_hedges:
					if h in v.hedgelist:

						v.hedgelist.remove(h)

			# Add new split hedges to the vertex hedgelist
			v.hedgelist += [h for h in new_hedges if h.origin == v]

		for v in update_vertices:
			print('v', v)

			l = len(v.hedgelist)

			if l < 2:
				raise DcelError(
					"Badly formed dcel: less than two hedges in vertex")
			elif l == 2:
				belong = siteBelong(new_site, close_site,
									v.hedgelist[0], v.hedgelist[1])
				for h in v.hedgelist:
					print(h)
					site = belong[h]
					#Fix here
					if lefton(h, site):
						tail = h
					else:
						head = h.twin
				tail.nexthedge = head
				head.prevhedge = tail
				#print('cur', tail)
				#print('next', tail.nexthedge)

			else:
				# Handle the new intersect vertices
				if v in intersect_vl:
					for h in v.hedgelist:
						print(h)

					v.sortthree(new_site, close_site, hb1)
				else:
    				#if mergevertex:
    					
					for h in v.hedgelist:
						if h in new_hedges:
							v.sortthree(new_site, close_site)
				
				pivothedge = v.hedgelist[0]
				print('pivothedge', pivothedge)
				p = pivothedge.origin.coord
				q = pivothedge.v1.coord
				compare_hedges = v.hedgelist[1:]

				#v.hedgelist[0].nexthedge = v.hedgelist[1].twin
				#print('cur',v.hedgelist[0])
				#print('next',v.hedgelist[0].nexthedge)
#
				#v.hedgelist[1].twin.prevhedge = v.hedgelist[0]

				#if len(v.hedgelist) == 2 or len(v.hedgelist) == 3:

				eps = 10e-3
				minAngle = float('inf')
				
				for h in compare_hedges:
					r = h.v1.coord
					print('ccw', ccw(p, q, r), pivothedge, h)
					if ccw(p, q, r) < 0 and ccw(p, q, r) < minAngle:
						left = h
						minAngle = ccw(p, q, r)
					
  			
					else:
						right = h 
						print('right', right)
				print('left', left)
				
				
				
				pivothedge.nexthedge = left.twin
				left.twin.prevhedge = pivothedge

				print('cur', pivothedge)
				print('next', pivothedge.nexthedge)

				p = right.origin.coord 
				q = right.v1.coord 
				r = left.v1.coord 
				
				if ccw(p, q, r) > 0 and np.abs(ccw(p, q, r)) > eps:
					right.twin.prevhedge = left
					left.nexthedge = right.twin
					print('cur', left)
					print('next', left.nexthedge)

					right.nexthedge = pivothedge.twin
					pivothedge.twin.prevhedge = right 
					print('cur', right)
					print('next', right.nexthedge)

				else:
					right.nexthedge = pivothedge.twin
					pivothedge.twin.prevhedge = right
					print('cur', right)
					print('next', right.nexthedge)


				
						
				#else:
				#	v.hedgelist.pop()

		# Step 4: Face assignment
		belong = siteBelong(new_site, close_site,
							newface_hedges[0], newface_hedges[1])

		for h in newface_hedges:


			site = belong[h]
			print('belong', site)

			if site == close_site:
				f = self.faces[site]
				f.hedges = []
				

			else:
				if ignoreface:
					continue
				# If this is a new site, create a new newface for that site
				f = Face(site=site)
				self.faces[f.site] = f
				# print('new', f.site)

			# We link the hedge to the new newface

			f.wedge = h
			f.hedges.append(h)
			f.wedge.newface = f
			f.vertices = []

			# And we traverse the boundary of the newface

			print(h)
			f.vertices.append(h.origin)
			linkheges = {h.v1: [h]}

			i = 0 
			
			while (not h.nexthedge is f.wedge) and i<10:
			
				h = h.nexthedge
				f.vertices.append(h.origin)
				f.hedges.append(h)

				if mergevertex:
					if h.v1 in linkheges:
						if h not in linkheges[h.v1]:
							linkheges[h.v1].append(h)
					else:
						linkheges[h.v1] = [h]


				print('next:', h)
				h.newface = f
				i += 1
			
			if i == 10:
				print('Equal 10')
				for key, value in linkheges.items():
					
					if len(value) == 2:

						head = value[0]
						tail = value[1]


						if head.prevhedge.twin != tail:
							head = value[1]
							tail = value[0]

						
						f.vertices.remove(key)
						f.vertices.remove(tail.origin)
						
						
						
						f.hedges.remove(head)
						f.hedges.remove(tail)

						tail= tail.prevhedge 
						print('head', head)
						print('tail', tail)

						newhedge1 = Hedge(tail.v1, head.origin)
						newhedge2 = Hedge(head.origin, tail.v1)
						newhedge1.twin = newhedge2
						newhedge2.twin = newhedge1
						
						newhedge1.nexthedge = head.nexthedge 
						head.nexthedge.prevhedge = newhedge1 
						print('cur', newhedge1)
						print('next', newhedge1.nexthedge)

						newhedge1.prevhedge = tail.prevhedge
						tail.prevhedge.nexthedge = newhedge1
						print('cur', newhedge1.prevhedge)
						print('next', newhedge1)
						
						newhedge1.newface = f
						f.hedges.append(newhedge1)

						
						f.vertices[f.vertices.index(tail.v1)].hedgelist.remove(tail.twin)
						f.vertices[f.vertices.index(tail.v1)].hedgelist.append(newhedge1.twin)
						print('remove',tail.twin)
						print('append', newhedge1.twin)


						f.vertices[f.vertices.index(head.origin)].hedgelist.remove(head)
						f.vertices[f.vertices.index(head.origin)].hedgelist.append(newhedge1)

						if f.wedge == head or f.wedge == tail:
							f.wedge = newhedge1
						
						
						break
			
				h = f.wedge
				print('cur', f.wedge)
				i = 0 
				while h.nexthedge != f.wedge and i <7:
					h = h.nexthedge 
					print('next',h)
					i += 1
				f.hedges = list(set(f.hedges))
				f.vertices = list(set(f.vertices))
						

    				

			
		# Return the hedge that is not the border

		return ignoreface, {v: h for v, h in intersect_edges.items() if not isborder(self.border, h) and h != intersected_hed}
	def findpoints(self, pl, onetoone=False):
		"""Given a list of points pl, returns a list of
		with the corresponding newface each point belongs to and
		None if it is outside the map.
		"""

		ans = []
		if onetoone:
			fl = self.faces[:]
			for p in pl:
				found = False
				for f in fl:
					if f.external:
						continue
					if f.isinside(p):
						fl.remove(f)
						found = True
						ans.append(f)
						break
				if not found:
					ans.append(None)

		else:
			for p in pl:
				found = False
				for f in self.faces:
					if f.external:
						continue
					if f.isinside(p):
						found = True
						ans.append(f)
						break
				if not found:
					ans.append(None)

		return ans

	def getFace(self, site):
		return self.faces[site]


def hsort(h1, h2):
	"""Sorts two half edges counterclockwise"""

	if h1.angle < h2.angle:
		return -1
	elif h1.angle > h2.angle:
		return 1
	else:
		return 0


def checkhedges(hl):
	"""Consistency check of a hedge list: nexthedge, prevhedge"""

	for h in hl:
		if h.nexthedge not in hl or h.prevhedge not in hl:
			raise DcelError("Problems with an orphan hedge..")


def area2(hedge, point):
	"""Determines the area of the triangle formed by a hedge and
	an external point"""

	pa = hedge.twin.origin
	pb = hedge.origin
	pc = point
	return (pb.x - pa.x)*(pc[1] - pa.y) - (pc[0] - pa.x)*(pb.y - pa.y)


def lefton(hedge, point):
	"""Determines if a point is to the left of a hedge"""
	#print(area2(hedge, point))

	return area2(hedge, point) >= 0


def hangle(dx, dy):
	"""Determines the angle with respect to the x axis of a segment
	of coordinates dx and dy
	"""

	l = m.sqrt(dx*dx + dy*dy)
	if dy > 0:
		return m.acos(dx/l)
	else:
		return 2*m.pi - m.acos(dx/l)


def split_hedge(v, hedge):
	v1 = hedge.v1
	v2 = hedge.origin

	# Edge with v is the tail
	htail1 = Hedge(v, v2)
	# Twin of htail1
	htail2 = Hedge(v2, v)
	# Edge with v is the origin
	horigin1 = Hedge(v1, v)
	# Twin of horigin
	horigin2 = Hedge(v, v1)
	return htail1, htail2, horigin1, horigin2


def minDistance(A, B, E):
	'''Function to return the minimum distance
 between a line segment AB and a point E'''

	# vector AB
	AB = [None, None]
	AB[0] = B[0] - A[0]
	AB[1] = B[1] - A[1]

	# vector BP
	BE = [None, None]
	BE[0] = E[0] - B[0]
	BE[1] = E[1] - B[1]

	# vector AP
	AE = [None, None]
	AE[0] = E[0] - A[0]
	AE[1] = E[1] - A[1]

	# Variables to store dot product

	# Calculating the dot product
	AB_BE = AB[0] * BE[0] + AB[1] * BE[1]
	AB_AE = AB[0] * AE[0] + AB[1] * AE[1]

	# Minimum distance from
	# point E to the line segment
	reqAns = 0

	# Case 1
	if (AB_BE > 0):

		# Finding the magnitude
		y = E[1] - B[1]
		x = E[0] - B[0]
		reqAns = sqrt(x * x + y * y)

	# Case 2
	elif (AB_AE < 0):
		y = E[1] - A[1]
		x = E[0] - A[0]
		reqAns = sqrt(x * x + y * y)

	# Case 3
	else:

		# Finding the perpendicular distance
		x1 = AB[0]
		y1 = AB[1]
		x2 = AE[0]
		y2 = AE[1]
		mod = sqrt(x1 * x1 + y1 * y1)
		reqAns = abs(x1 * y2 - y1 * x2) / mod

	return reqAns

def siteClose(site1, site2, hedge):

	A = hedge.vertices[0].coord
	B = hedge.vertices[1].coord
	if minDistance(A, B, site1) < minDistance(A, B, site2):
		return site1 
	else:
		return site2

def siteBelong(site1, site2, hedge1, hedge2):

	hedges = [hedge1, hedge2]
	belong = {}
	distance = {}
	alternative = {}

	for hedge in hedges:
		A = hedge.vertices[0].coord
		B = hedge.vertices[1].coord
		if minDistance(A, B, site1) < minDistance(A, B, site2):
			belong[hedge] = site1
			distance[hedge] = minDistance(A, B, site1)
			alternative[hedge] = site2

		else:
			belong[hedge] = site2
			distance[hedge] = minDistance(A, B, site2)
			alternative[hedge] = site1

	# If both hedges are close to the same site, the hedge with closer distance to the site will belong to that site
	if belong[hedges[0]] == belong[hedges[1]]:
		if distance[hedges[0]] < distance[hedges[1]]:
			belong[hedge2] = alternative[hedges[1]]
		else:
			belong[hedge1] = alternative[hedges[0]]

	return belong


def todelete_hedge(hedge1, hedge2, point):

	if minDistance(hedge1.v1.coord, hedge1.origin.coord, point) < minDistance(hedge2.v1.coord, hedge2.origin.coord, point):
		return hedge1

	else:
		return hedge2


def isborder(border, hedge):

	if hedge.origin.coord[0] == hedge.v1.coord[0] and hedge.origin.coord[0] in border:
		return True

	if hedge.origin.coord[1] == hedge.v1.coord[1] and hedge.origin.coord[1] in border:
		return True

	else:
		return False

def isOnLine(pt, hedge):

	pt1 = hedge.v1.coord 
	pt2 = hedge.origin.coord 
	pt3 = pt 

	x1, x2, x3 = pt1[0], pt2[0], pt3[0]

	y1, y2, y3 = pt1[1], pt2[1], pt3[1]

	if x1 == x2:
		if x1 == x3:
			pt3_on = True
		else:
			pt3_on = False

	else:
		slope = (y2 - y1) / (x2 - x1)

		pt3_on = (y3 - y1) == slope * (x3 - x1)

	pt3_between = (min(x1, x2) <= x3 <= max(x1, x2)) and (
		min(y1, y2) <= y3 <= max(y1, y2))

	on_and_between = pt3_on and pt3_between

	return on_and_between

#the representation of a vector


class vec:
	def __init__(self, _x, _y):
		self.x = _x
		self.y = _y

#multiply vector by scalar


def scale(v, s):
	return vec(s*v.x, s*v.y)

#point to vector


def toVec(p1, p2):
	return vec(p2[0]-p1[0], p2[1]-p1[1])

#translate a point

#return the dot product


def dot(v1, v2):
	return v1.x*v2.x + v1.y*v2.y

def norm(v):
	return (v.x*v.x + v.y*v.y)**0.5


def angle(v1, v2):
	v1_norm = norm(v1)
	v2_norm = norm(v2)
	cos = dot(v1, v2)/(v1_norm*v2_norm)
	return np.arccos(cos)

#cross product

def cross(v1, v2):
	return v1.x*v2.y - v2.x*v1.y

#return true if the point r is on the left side of pq

def ccw(p, q, r):

	return cross(toVec(p, q), toVec(p, r))
	
def is_right(h1, h2):
	right = False
	v1 = toVec(h1.v1.coord, h1.origin.coord)
	v2 = toVec(h2.v1.coord, h2.origin.coord)
	dotprod = dot(v1, v2)
	print(dotprod)
	if dotprod > 0:
		right = True 
	return right
			
			
	
	
	
		

