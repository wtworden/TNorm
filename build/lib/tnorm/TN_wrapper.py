#-*-python-*-
from __future__ import print_function

import regina
import snappy
import sys
import time


from tnorm.euler import *
from tnorm.boundary import *
from tnorm.homology import *
from tnorm.regina_helpers import *
from tnorm.norm_ball import *
from tnorm.sage_types import *
from tnorm.matrices import *
from tnorm.cached_property import *

import tnorm

#preparser(False)


class TN_wrapper():
	"""
	A class for working with Thurston Norm on H2(M,\partial M).

	sage: W = TN_wrapper(SnapPeaTri)

	SnapPeaTri can be any of the following:
	(1) A snappy Manifold
	(2) a regina.SnapPeaTriangulation triangulation
	(3) Any string argument that snappy.Manifold() can take, e.g.,

	sage: W = TN_wrapper('L6a5')
	sage: W = TN_wrapper('/path/to/snappy/file.tri')
	sage: W = TN_wrapper(string)       # where string is the triangulation encoded as a string (i.e., contents of a tri file)

	"""
	def __init__(self,SnapPeaTri,qtons=None,tracker=False,quiet=False,allows_non_admissible=False):
		if isinstance(SnapPeaTri,str):
			self.manifold = snappy.Manifold(SnapPeaTri)
			self.triangulation = regina.SnapPeaTriangulation(self.manifold._to_string())
		elif isinstance(SnapPeaTri, regina.engine.SnapPeaTriangulation):
			self.triangulation = SnapPeaTri
			self.manifold = snappy.Manifold(self.triangulation.snapPea())
		elif isinstance(SnapPeaTri, snappy.Manifold):
			self.manifold = SnapPeaTri
			self.triangulation = regina.SnapPeaTriangulation(SnapPeaTri._to_string())
		self._angle_structure = solve_lin_gluingEq(self.triangulation)
		self._qtons = qtons
		self._tkr = False
		self.allows_non_admissible = allows_non_admissible
		tnorm.QUIET = quiet

		# qtons memoization caches-----
		self._euler_char = {}
		self._boundary_slopes = {}
		self._map_to_H1bdy = {}
		self._map_to_ball = {}
		self._map_to_H2 = {}
		self._num_boundary_comps_H1 = {}
		self._num_boundary_comps = {}
		self._over_face = {}
		self._is_norm_minimizing = {}

		if not tnorm.QUIET:			
			print('Enumerating quad transversely oriented normal surfaces (qtons)... ', end='')
		if tracker:
			self._tkr = regina.ProgressTracker()
		self._qtons = self.qtons()
		for i in range(self._qtons.size()):
			self._qtons.surface(i).setName(str(i))
		if not tnorm.QUIET:
			print('Done.\n')

	def qtons(self):
		"""
		Enumerate vertex oriented quad normal surfaces for the ideal triangulation self.triangulation()
		This is cached in self.qtons, so calls to this function after the initial call (when TN_wrapper() is 
		instantiated) will be fast. 

		Returns a regina.NormalSurfaces object.
		"""
		if self._qtons == None:
			if self._tkr:
				if self.allows_non_admissible:
					return regina.NormalSurfaces.enumerate(self.triangulation,regina.NS_ORIENTED_QUAD,regina.NS_IMMERSED_SINGULAR,regina.NS_ALG_DEFAULT,self._tkr)
				else:
					return regina.NormalSurfaces.enumerate(self.triangulation,regina.NS_ORIENTED_QUAD,regina.NS_VERTEX,regina.NS_ALG_DEFAULT,self._tkr)
			else:
				if self.allows_non_admissible:
					return regina.NormalSurfaces.enumerate(self.triangulation,regina.NS_ORIENTED_QUAD,regina.NS_IMMERSED_SINGULAR,regina.NS_ALG_DEFAULT)
				else:
					return regina.NormalSurfaces.enumerate(self.triangulation,regina.NS_ORIENTED_QUAD,regina.NS_VERTEX,regina.NS_ALG_DEFAULT)
		else:
			return self._qtons


	def show_progress(self):
		print(self._tkr.percent())


	def euler_char(self,sn_surface):
		"""
		Return the Euler characteristic of the spun normal surface. 
		Argument can be a surface (Regina oriented spun normal surface) or the index of a surface.

		sage: W.euler_char(3)    # return Euler characteristic of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		-1

		"""
		try:
			ind = int(sn_surface)
			s = self.qtons().surface(ind)
		except TypeError:
			ind = int(sn_surface.name())
			s = sn_surface
			
		if ind in self._euler_char:
			return self._euler_char[ind]
		else:
			ec = euler_char_(s,self._angle_structure)
			self._euler_char[ind] = ec
			return ec


	def boundary_slopes(self,sn_surface):
		"""
		Return the boundary slopes of the spun normal surface. This will be different than the return value
		of Regina's method boundaryIntersections(), because Regina returns the intersection numbers,
		not the slope (i.e., Regina's method will be the negative inverse of these slopes). 
		Argument can be a surface (Regina oriented spun normal surface) or the index of a surface.

		sage: W.boundary_slopes(3)    # return boundary slopes of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		[(-1, 0), (0, 1), (-1, 0)]

		"""
		try:
			ind = int(sn_surface)
			s = self.qtons().surface(ind)
		except TypeError:
			ind = int(sn_surface.name())
			s = sn_surface

		if ind in self._boundary_slopes:
			return self._boundary_slopes[ind]
		else:
			bs = boundary_slopes_(s, self.manifold)
			self._boundary_slopes[ind] = bs
			return bs

	def map_to_H1bdy(self,sn_surface):
		"""
		Return the image of the given surface in H1(\partial M). Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		sage: W._map_to_H1bdy(3)    # return image in H1(\partial M) of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		[(-1, 0), (0, 1), (-1, 0)]

		"""
		try:
			ind = int(sn_surface)
			s = self.qtons().surface(ind)
		except TypeError:
			ind = int(sn_surface.name())
			s = sn_surface

		if ind in self._map_to_H1bdy:
			return self._map_to_H1bdy[ind]
		else:
			h1 = qtons_to_H1bdy(s, self.manifold)
			self._map_to_H1bdy[ind] = h1
			return h1
	
	def map_to_H2(self,sn_surface):
		"""
		Return the image of the given surface in H2(M, \partial M). Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		sage: W.mapToH2(3)    # return the image in H2(M,\partial M) of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		(0,1,0)

		"""
		try:
			ind = int(sn_surface)
			s = self.qtons().surface(ind)
		except TypeError:
			ind = int(sn_surface.name())
			s = sn_surface

		if ind in self._map_to_H2:
			return self._map_to_H2[ind]
		else:
			h2 = homology_map(s, self.manifold)
			self._map_to_H2[ind] = h2
			return h2

	def over_face(self, sn_surface):

		try:
			ind = int(sn_surface)
		except TypeError:
			ind = int(sn_surface.name())

		if ind in self._over_face:
			return self._over_face[ind]
		else:
			v = self.map_to_H2(ind)
			f = over_face_(v,self.TNormBall.polyhedron)
			self._over_face[ind] = f
			return f

	def is_norm_minimizing(self, sn_surface):

		try:
			ind = int(sn_surface)
		except TypeError:
			ind = int(sn_surface.name())

		if ind in self._is_norm_minimizing:
			return self._is_norm_minimizing[ind]
		else:
			v = self.map_to_H2(ind)
			if not self.TNormBall.polyhedron.interior_contains(v):
				self._is_norm_minimizing[ind] = True
				return True
			else:
				self._is_norm_minimizing[ind] = False
				return False

	def map_to_ball(self,sn_surface):
		"""
		Return the image of the given surface in H2(M, \partial M), then divide by Euler characteristic. Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		"""
		try:
			ind = int(sn_surface)
		except TypeError:
			ind = int(sn_surface.name())

		if ind in self._map_to_ball:
			return self._map_to_ball[ind]
		else:
			image_in_H2 = self.map_to_H2(ind)
			if self.euler_char(ind) < 0:
				mtb = tuple([Rational(i)/(-self.euler_char(ind)) for i in image_in_H2])
				self._map_to_ball[ind] = mtb
				return mtb
			else:
				self._map_to_ball[ind] = Infinity
				return Infinity



	def genus(self,sn_surface):
		"""
		Return the genus. Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		sage: W.genus(3)    # return genus of the surface with index 3, i.e, W.qtons().surface(3).
		1
		"""
		return (2-self.euler_char(sn_surface)-self.num_boundary_comps(sn_surface))/2
	
	def num_boundary_comps(self,sn_surface):
		"""
		Return the number of boundary components. Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		sage: W.num_boundary_comps(3)    # return num of boundary components of surface with index 3, i.e, W.qtons().surface(3).
		2

		"""	
		try:
			ind = int(sn_surface)
			s = self.qtons().surface(ind)
		except TypeError:
			ind = int(sn_surface.name())
			s = sn_surface

		if ind in self._num_boundary_comps:
			return self._num_boundary_comps[ind]
		else:
			nbc = num_boundary_comps_(s)
			self._num_boundary_comps[ind] = nbc
			return nbc


	def num_boundary_comps_H1(self,sn_surface):
		"""
		Return the number of boundary components after mapping to H1(bdy m). This will be the number
		of boundary components after cancelling oppositely oriented pairs (i.e., a 3-puntured sphere
		will become a once-punctured torus if two of the punctures are on the cusps and boundary curves
		with opposite orientations). Argument can be a surface (regina oriented spun normal surface) or 
		the index of a surface.

		sage: W.numBoundaryComps(3)    # return num of boundary components of surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		2

		"""	
		try:
			ind = int(sn_surface)
			s = self.qtons().surface(ind)
		except TypeError:
			ind = int(sn_surface.name())
			s = sn_surface

		if ind in self._num_boundary_comps_H1:
			return self._num_boundary_comps_H1[ind]
		else:
			nbc = 0
			b = self.map_to_H1bdy(ind)
			for slope in b:
				nbc += gcd(slope[0],slope[1])
			self._num_boundary_comps_H1[ind] = nbc
			return nbc

		
	def _image_in_H2M(self):
		image = []
		for i in range(self.qtons().size()):
			pt = self.map_to_H2(i)
			if not pt.is_zero():
				ec = self.euler_char(i)
				if (pt,ec) not in [tup[1:] for tup in image]:
					image.append((i,pt,ec))
		return image

	@cached_property
	def norm_ball_points(self):
		"""
		Return the points in H2(M,bdy M) which are images of (normalized) vertices of the oriented spun
		normal surface projective solution space with negative Euler characteristic. The convex hull of 
		these points is the norm ball.
		"""
		image = dict()
		image_on_ball = []
		image_w_euler = self._image_in_H2M()
		image_w_euler.sort(key=lambda x: x[0])
		image_w_euler.reverse()
		for (i,pt,ec) in image_w_euler:
			if ec < 0:
				image_on_ball.append(tuple([coord/(-ec) for coord in pt]))
				image[tuple([coord/(-ec) for coord in pt])]=i
		return image

	@cached_property
	def TNormBall(self):
		"""
		Return the norm ball as a sage Polyhedron. Some things you can do:

		sage: B=W.TNormBall
		sage: B.faces(2)        # list the 2-dimensional faces
		sage: B.faces(0)        # list the 0-dimensional faces (vertices)

		sage: F = B.faces(2)
		sage: F[0].vertices()   #return vertices of the 0 face

		sage: B.plot()          # plot of B. I imagine this only works if B is 1,2, or 3 dimensional.

		The Polyhedron class has a lot of methods, and I have not explored many of them. There is 
		probably a lot more that would be useful. Use tab completion to explore!
		"""
		if not tnorm.QUIET:
			print('Computing Thurston norm unit ball... ', end='')
		pts_dict = self.norm_ball_points
		poly = Polyhedron(vertices=Matrix(pts_dict.keys()))
		verts = tuple([(pts_dict[tuple(v)],vector(v)) for v in poly.vertices_list()])
		vertices = [Vertex(i,verts[i][0],verts[i][1],self.num_boundary_comps_H1(verts[i][0]),self.euler_char(verts[i][0]),self.boundary_slopes(verts[i][0]),self.map_to_H1bdy(verts[i][0])) for i in range(len(verts))]
		ball = NormBall(vertices, poly)
		if not tnorm.QUIET:
			print('Done.\n')
		return ball


	@cached_property
	def qtons_info(self):
		if not tnorm.QUIET:
			print('Analyzing quad transversely oriented normal surfaces... ', end='')
		qtons_info_dict = dict([(i,{}) for i in range(self.qtons().size())])
		for i in range(self.qtons().size()):
			qtons_info_dict[i]['image_in_H2'] = self.map_to_H2(i)
			qtons_info_dict[i]['euler_char'] = self.euler_char(i)
			qtons_info_dict[i]['qtons_slopes'] = self.boundary_slopes(i)
			qtons_info_dict[i]['H1_boundary_slopes'] = self.map_to_H1bdy(i)
			qtons_info_dict[i]['num_boundary_comps'] = self.num_boundary_comps(i)
			qtons_info_dict[i]['genus'] = (2-self.euler_char(i)-self.num_boundary_comps(i))/2
			qtons_info_dict[i]['is_norm_minimizing'] = self.is_norm_minimizing(i)
			qtons_info_dict[i]['over_face'] = self.over_face(i)
		if not tnorm.QUIET:
			print('Done.\n')
		return qtons_info_dict

def check_subfaces(ray,f_poly):
	if f_poly.dim() == 0:
		return f_poly
	for f in f_poly.faces(f_poly.dim()-1):
		if not f.as_polyhedron().intersection(ray).is_empty():
			return check_subfaces(ray,f.as_polyhedron())
		else:
			continue
	return f_poly
			

def over_face_(v, P):
	if v.is_zero():
		return None
	else:
		ray = Polyhedron(rays=[v])
		f_poly = check_subfaces(ray, P)
		for i in range(P.dim()):
			for f in P.faces(i):
				if sorted([vert.vector() for vert in f.vertices()]) == sorted([vert.vector() for vert in f_poly.vertices()]):
					return f





















