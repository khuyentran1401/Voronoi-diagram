class Point:
	def __init__(self, x, y):
		self.x = x
		self.y = y

	def __str__(self): 
		return str((self.x, self.y))

	def __repr__(self):
		return (self.x, self.y)


class Line:
	def __init__(self, m, b):
		self.m = m
		self.b = b


def midpoint(p1, p2):

	m_p = ((p1[0] + p2[0])/2, (p1[1] + p2[1])/2)

	return m_p


def find_b(p, m):
	return p[1] - m * p[0]


def slope(p1, p2):

	# If the slope is vertical
	if p2[0] == p1[0]:
		return 'undefined'
	return (p2[1] - p1[1])/(p2[0] - p1[0])


def slope_perpendicular_bisector(p1, p2):
		
	# If the line is horizontal
	if p2[0] == p1[0]:
		return 0
	
	if p2[1] == p1[1]:
		return 'undefined'
		
	return -((p2[0] - p1[0])/(p2[1] - p1[1]))


def perpendicular_bisector(p1, p2, xmin, xmax, ymin, ymax):


	m = slope_perpendicular_bisector(p1, p2)
	
	m_p = midpoint(p1, p2)

	if m == 'undefined':
		p_r1 = (m_p[0], ymin-2 )
		p_r2 = (m_p[0], ymax+2 )
			
	elif m == 0:
		p_r1 = (xmin-2, m_p[1])
		p_r2 = (xmax+2, m_p[1])
	
	else:

		b = find_b(m_p, m)

		x_small = xmin - 2
		x_large = xmax + 2

		# Find points that the bisector intersect the rectangle
		# Add more cases
		y_r1 = m*(x_small - 1) + b
		p_r1 = (x_small - 1, y_r1)
		y_r2 = m*(x_large+1) + b
		p_r2 = (x_large + 1, y_r2)

	return p_r1, p_r2
	
def intersection(p1, q1, p2, q2, shift=None):
		

	m1 = slope(p1, q1)
	if shift:
		m1 -= shift

	m2 = slope(p2, q2)

	if m1 != 'undefined':
		b1 = find_b(p1, m1)
	if m2 != 'undefined':
		b2 = find_b(p2, m2)

	# If the line is vertical, we take the x of the point to be b

	if m1 == 'undefined' and m2 == 0:
		x = p1[0]
		y = p2[1]
		
	elif m2 == 'undefined' and m1 == 0:
		x = p2[0]
		y = p1[1]

	elif m1 == 'undefined':
		x = p1[0]
		y = m2*x+b2

	elif m2 == 'undefined':
		x = p2[0]
		y = m1*x+b1
		

	else:
		x = (b2-b1)/(m1-m2)
	
		y = m1*x+b1

	return (x, y)


if __name__ == '__main__':
	p1 = (7, -5)
	p2 = (-7, 5)

	print(slope(p1, p2))
