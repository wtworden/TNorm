#-*-python-*-
from __future__ import print_function

import sys
try:
    sys.path.remove('/Applications/Regina.app/Contents/MacOS/python')
except ValueError:
    pass
import regina

import snappy
import time
import itertools


#import imp
#regina = imp.load_source('regina', '/Applications/SageMath/local/lib/python2.7/site-packages/sageRegina-5.1.5-py2.7.egg-info')

from tnorm.kernel.simplicial import get_face_map_to_C2, get_quad_map_to_C2, H2_as_subspace_of_C2, qtons_image_in_C2
from tnorm.kernel.boundary import bdy_slopes_unoriented_
from tnorm.kernel.euler import solve_lin_gluingEq, euler_char_
from tnorm.kernel.homology import qtons_to_H1bdy, homology_map
from tnorm.kernel.matrices import intersection_mat, oriented_quads_mat
from tnorm.kernel.regina_helpers import regina_to_sage_int
from tnorm.norm_ball import *
from tnorm.sage_types import *
from tnorm.utilities import *
from tnorm.x3d_to_html import *
import tnorm.constants



#from sympy import Matrix
#from sympy import Rational as QQ


#preparser(False)


class TN_wrapper():
	"""
	A class for working with Thurston Norm on H2(M,\partial M).

	sage: W = TN_wrapper(SnapPeaTri)

	SnapPeaTri can be any of the following:
	(1) A snappy Manifold or Triangulation
	(2) a regina.SnapPeaTriangulation triangulation
	(3) Any string argument that snappy.Manifold() can take, e.g.,

	sage: W = TN_wrapper('L6a5')
	sage: W = TN_wrapper('/path/to/snappy/file.tri')
	sage: W = TN_wrapper(string)       # where string is the triangulation encoded as a string (i.e., contents of a tri file)
	
	"""
	def __init__(self, manifold, qtons=None, quiet=False, tracker=False, allows_non_admissible=False, bdy_H1_basis='natural', force_simplicial_homology=False):
		
		if isinstance(manifold,str):
			self.manifold = snappy.Manifold(manifold)
		elif isinstance(manifold, regina.engine.SnapPeaTriangulation):
			self.manifold = snappy.Manifold(manifold.snapPea())
		elif isinstance(manifold, regina.engine.Triangulation3):
			self.manifold = snappy.Manifold(manifold.snapPea())
		elif isinstance(manifold, snappy.Manifold):
			self.manifold = manifold
		elif isinstance(manifold, snappy.Triangulation):
			self.manifold = manifold

		for c in self.manifold.cusp_info():
			if c.is_complete == False:
				self.manifold = self.manifold.filled_triangulation()
				break
		
		self._QUIET = quiet
		self._force_simplicial_homology = force_simplicial_homology
		self.num_cusps = self.manifold.num_cusps()
		if self.num_cusps != 0:
			try: 
				L = self.manifold.link()
				self.knows_link_complement = True
			except ValueError:
				self.knows_link_complement = False
			if self.knows_link_complement and bdy_H1_basis == 'natural':
				self.bdy_H1_basis = 'natural'
			else:
				self.manifold.set_peripheral_curves('shortest')
				self.bdy_H1_basis = 'shortest'

			self.triangulation = regina.SnapPeaTriangulation(self.manifold._to_string())
			self._angle_structure = solve_lin_gluingEq(self.triangulation)
			self._intersection_matrices = [intersection_mat(self.manifold, self.triangulation, cusp) for cusp in range(self.triangulation.countCusps())]
			self.manifold_is_closed = False
		else:
			self.triangulation = regina.Triangulation3(self.manifold._to_string())
			self.triangulation.intelligentSimplify()
			self.manifold_is_closed = True
			self.bdy_H1_basis = None
		if not self.triangulation.isOriented():
			self.triangulation.orient()
		self._qtons = qtons
		self._tkr = False
		self.allows_non_admissible = allows_non_admissible
		self.is_fibered = 'unknown'
		self.betti_number = self.triangulation.homologyH1().rank()

		# qtons memoization caches-----
		self._euler_char = {}
		self._map_to_ball = {}
		self._map_to_H2 = {}
		self._num_boundary_comps = {}
		self._over_facet = {}
		self._is_norm_minimizing = {}
		self._is_admissible = {}
		self._qtons_image_in_C2 = {}

		if self.manifold.num_cusps() > 0:
			self._map_to_H1bdy = {}
			self._regina_bdy_slopes = {}

		if not self._QUIET:		
			print('Enumerating quad transversely oriented normal surfaces (qtons)... ', end='')
			try:
				sys.stdout.flush()
			except AttributeError:
				pass

		if tracker:
			self._tkr = regina.ProgressTracker()

		# compute transversely oriented normal surfaces
		self._qtons = self.qtons()

		# name the surfaces by their index in the NormalSurfaces list
		for i in range(self._qtons.size()):
			self._qtons.surface(i).setName(str(i))

		if not self._QUIET:
			print('Done.')
			try:
				sys.stdout.flush()
			except AttributeError:
				pass

		if self.betti_number > self.num_cusps:
			self.has_internal_homology = True
		else:
			self.has_internal_homology = False

		# if the manifold has internal homology or force_simplicial_homology==True then we need to use simplicial homology.
		if self.has_internal_homology or self._force_simplicial_homology:
			
			self.uses_simplicial_homology = True

			if not self._QUIET:
				print('computing simplicial homology...',end='')
				try:
					sys.stdout.flush()
				except AttributeError:
					pass

			self._face_map_to_C2 = get_face_map_to_C2(self.triangulation)
			self._quad_map_to_C2 = get_quad_map_to_C2(self.triangulation, self._face_map_to_C2)
			H2_basis_in_C2, P, qtons_image = H2_as_subspace_of_C2(self, self._face_map_to_C2, self._quad_map_to_C2)
			self._project_to_im_del3 = P
			self._qtons_image_in_C2 = {i:qtons_image[i] for i in range(len(qtons_image))}
			assert len(H2_basis_in_C2) == self.betti_number, self.manifold.name()+', force_simplicial_homology={}'.format(self.force_simplicial_homology)
				# raise an error if the qtons do not generate H2. This should only happen for a knot in a 
				# rational homology sphere, with force_simplicial_homology=True.
			I = Matrix.identity(self.betti_number)
			B = Matrix(H2_basis_in_C2).transpose()
			A = B.solve_left(I)
			self._map_H2_to_standard_basis = A

			if not self._QUIET:
				print('Done.')
				try:
					sys.stdout.flush()
				except AttributeError:
					pass
		else:
			self.uses_simplicial_homology = False

	# need to make this work for qtons without index. I.e., a linear combo of vertex normal surfaces.
	def _simplicial_map_to_H2(self,qtons):
		if self.uses_simplicial_homology:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
			if ind in self._map_to_H2:
				return self._map_to_H2[ind]
			else:
				s_H2 = self._map_H2_to_standard_basis*(self._qtons_image_in_C2[ind]-self._project_to_im_del3*self._qtons_image_in_C2[ind])
				self._map_to_H2[ind] = s_H2
				return s_H2
		else:
			return None

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
		try:
			sys.stdout.flush()
		except AttributeError:
			pass


	def euler_char(self,qtons):
		"""
		Return the Euler characteristic of the normal surface. 
		Argument can be a surface (Regina oriented spun normal surface) or the index of a surface.

		sage: W.euler_char(3)    # return Euler characteristic of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		-1

		"""
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())
			
		if ind in self._euler_char:
			return self._euler_char[ind]
		else:
			s = self.qtons().surface(ind)
			if not self.manifold_is_closed:
				ec = euler_char_(s,self._angle_structure)
			else:
				ec = regina_to_sage_int(s.eulerChar())
			self._euler_char[ind] = ec
			return ec

	def boundary_slopes(self,qtons,quad_mat=None):
		"""
		Return the boundary slope of the given surface in H1(\partial M), with respect to the basis W.bdy_H1_basis. 
		Argument can be a surface (regina quad transverely oriented normal surface) or the index of a surface.

		Note that this slope may not be the same as W.regina_boundary_slopes(surface), as Regina computes slopes based on
		an orientation that depends on the direction that the spun normal surface spins into the cusp. Our orientation on 
		boundary slopes is instead based on the transverse orientation of the surface.

		sage: W._map_to_H1bdy(3)    # return image in H1(\partial M) of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		[(-1, 0), (0, 1), (-1, 0)]
		"""
		if self.manifold_is_closed:
			return []
		else:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if ind in self._map_to_H1bdy:
				return self._map_to_H1bdy[ind]
			else:
				s = self.qtons().surface(ind)
				h1 = qtons_to_H1bdy(s, self,quad_mat)
				self._map_to_H1bdy[ind] = h1
				return h1

	def regina_bdy_slopes(self,qtons,quad_mat=None):
		""" Return boundary slopes that result if we ignore transverse orientation, and cancel slopes that spin into the cusp
		in opposite directions. This may differ from slopes as computed by Regina by a sign, but should otherwise be the same.
		"""

		if self.manifold_is_closed:
			return []
		else:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if ind in self._regina_bdy_slopes:
				return self._regina_bdy_slopes[ind]
			else:
				s = self.qtons().surface(ind)
				rbs = bdy_slopes_unoriented_(s, self,quad_mat)
				self._regina_bdy_slopes[ind] = rbs
				return rbs

	def map_to_H2(self,qtons):
		"""
		Return the image of the given surface in H2(M, \partial M). Argument can be a surface (regina quad transversely oriented
		normal surface) or the index of a surface.

		sage: W.mapToH2(3)    # return the image in H2(M,\partial M) of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		(0,1,0)

		"""
		if not self.uses_simplicial_homology:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if ind in self._map_to_H2:
				return self._map_to_H2[ind]
			else:
				s = self.qtons().surface(ind)
				slopes = self.boundary_slopes(s)
				s_H2 = vector([slopes[i][1] for i in range(len(slopes))])
				self._map_to_H2[ind] = s_H2
				return s_H2
		else:
			return self._simplicial_map_to_H2(qtons)

	def simplicial_class(self,qtons):
		""" Return the vector in the dimension 2 chain group C2 corresponding to the qtons surface. The i^th coordinate
			of this vector corresponds to the number of triangles of index i, as indexed by Regina, with orientation positive
			if the transverse orientation of the surface points toward the tetrahedron triangle(i).front(). 
		"""
		if not self.uses_simplicial_homology:
			return None
		else:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if ind in self._qtons_image_in_C2:
				c = self._qtons_image_in_C2[ind]
			else:
				c = qtons_image_in_C2(self,self._quad_map_to_C2)
				self._qtons_image_in_C2[ind] = c
			return c

	def over_facet(self, qtons, as_string=False):
		"""Returns the facet of the norm ball that the qtons surface is above. E.g., <0 1 5> corresponds to the facet 
		with vertices 0, 1, and 5. 
		"""
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())

		if ind in self._over_facet:
			f = self._over_facet[ind]
		else:
			v = self.map_to_H2(ind)
			f = over_facet_(v,self.norm_ball.polyhedron)
			self._over_facet[ind] = f
		if as_string:
			return None if f==None else '<{}>'.format(' '.join([str(v.index()) for v in f.vertices()]))
		else:
			return f

	def is_norm_minimizing(self, qtons):
		"""Return True if the given qtons surface is Thurston norm minimizing.
		"""
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())

		if ind in self._is_norm_minimizing:
			return self._is_norm_minimizing[ind]
		else:
			v = self.map_to_H2(ind)
			if not self.norm_ball.polyhedron.interior_contains(v):
				self._is_norm_minimizing[ind] = True
				return True
			else:
				self._is_norm_minimizing[ind] = False
				return False

#	def locally_compatible(self,qtons1,qtons2):
#		pass

	def is_admissible(self, qtons):
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())

		if ind in self._is_admissible:
			return self._is_admissible[ind]
		else:
			s = self.qtons().surface(ind)
			mat = oriented_quads_mat(s)
			for row in mat:
				nonzero = [i for i in range(0,6,2) if ( row[i]!=0 or row[i+1]!=0 )]
				if len(nonzero) > 1:
					self._is_admissible[ind] = False
					return False
			self._is_admissible[ind] = True
			return True

	def map_to_ball(self, qtons):
		"""
		Return the image of the given surface in H2(M, \partial M), then divide by Euler characteristic. Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		"""
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())

		if ind in self._map_to_ball:
			return self._map_to_ball[ind]
		else:
			image_in_H2 = self.map_to_H2(ind)
			if self.euler_char(ind) < 0:
				mtb = tuple([QQ(i)/(-self.euler_char(ind)) for i in image_in_H2])
				self._map_to_ball[ind] = mtb
				return mtb
			else:
				self._map_to_ball[ind] = float('inf')
				return float('inf')

	def genus(self,qtons):
		"""
		Return the genus. Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		sage: W.genus(3)    # return genus of the surface with index 3, i.e, W.qtons().surface(3).
		1
		"""
		return (2-self.euler_char(qtons)-self.num_boundary_comps(qtons))/2

	def num_boundary_comps(self,qtons):
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
			ind = int(qtons)
			s = self.qtons().surface(ind)
		except TypeError:
			ind = int(qtons.name())
			s = qtons

		if ind in self._num_boundary_comps:
			return self._num_boundary_comps[ind]
		else:
			nbc = 0
			b = self.boundary_slopes(ind)
			for slope in b:
				nbc += gcd(slope[0],slope[1])
			self._num_boundary_comps[ind] = nbc
			return nbc

		
	def _image_in_H2M(self):
		image = []
		for i in range(self.qtons().size()):
			pt = self.map_to_H2(i)
			if not pt.is_zero():
				ec = self.euler_char(i)
				#if (pt,ec) not in [tup[1:] for tup in image]:
				image.append((i,pt,ec))
		return image

	@cached_property
	def _norm_ball_points(self):
		"""
		Return the points in H2(M,bdy M) which are images of (normalized) vertices of the oriented spun
		normal surface projective solution space with negative Euler characteristic. The convex hull of 
		these points is the norm ball.
		"""
		points0 = dict()
		points = dict()
		rays = dict()  # if the Euler char of s is 0 we can't normalize. This means the norm ball is
		               # infinite in the direction of map_to_H2(s), so it contains the ray (0,map_to_H2(s)).
		image_w_euler = self._image_in_H2M()


		for (i,pt,ec) in image_w_euler:
			if ec < 0:
				normalized = tuple([coord/(-ec) for coord in pt]) # divide by euler char to normalize

#				if not self.uses_simplicial_homology:
				if normalized not in points0:
					points0[normalized] = []
					points0[normalized].append((i,pt,ec))
				else:
					points0[normalized].append((i,pt,ec))

			elif ec == 0:
				rays[tuple(pt)] = i

		for nzd in points0:

			min_ec = min(-rep[2] for rep in points0[nzd])
			reps = [rep[0] for rep in points0[nzd] if -rep[2]==min_ec] # restrict to norm minimizing reps

			if self._force_simplicial_homology:
				# order reps by L^1 norm
				reps = sorted(reps,key=lambda x: sum(self.simplicial_class(x).apply_map(abs)))
				points[nzd] = reps[0] # get the index of the rep with minimal L^1 norm
			else:
				# order according to how well self.boundary_slopes(i) agrees with self.regina_bdy_slopes(i)
				if self.num_boundary_comps != 0:
					def sort_fn(rep):
						bs = self.boundary_slopes(rep)
						rbs = self.regina_bdy_slopes(rep)
						if bs == rbs:
							return 0
						elif all([rbs[i] in [bs[i],tuple([-bs[i][j] for j in range(len(bs[i]))])] for i in range(len(rbs))]):
							return 1
						else:
							return 2
					reps = sorted(reps,key=lambda x: sort_fn(x))
				points[nzd] = reps[0]

		return points, rays

	@cached_property
	def norm_ball(self):
		"""
		Return the norm ball as a sage Polyhedron. Some things you can do:

		sage: B=W.norm_ball
		sage: B.facets(2)        # list the 2-dimensional facets
		sage: B.facets(0)        # list the 0-dimensional facets (vertices)

		sage: F = B.facets(2)
		sage: F[0].vertices()   #return vertices of the 0^th facet

		sage: B.plot()          # plot of B. I imagine this only works if B is 1,2, or 3 dimensional.

		The Polyhedron class has a lot of methods, and I have not explored many of them. There is 
		probably a lot more that would be useful. Use tab completion to explore!
		"""
		if not self._QUIET:
			print('Computing Thurston norm unit ball... ', end='')
			try:
				sys.stdout.flush()
			except AttributeError:
				pass
		pts_dict, rays_dict = self._norm_ball_points

		### If M is not hyperbolic then rays_dict should be non-empty, and the norm ball should be non-compact. 
		### The below attempts to compute the non-compact norm ball, but it may be wrong. This is because (it seems) 
		### some surfaces may not be realized as spun normal surfaces.
		if len(rays_dict) != 0:
			polyhedron = Polyhedron(vertices=Matrix(pts_dict.keys()), rays=Matrix(rays_dict.keys()), base_ring=QQ, backend='cdd')
			rays = tuple([(rays_dict[tuple(i*vector(v))],i*vector(v)) for v in polyhedron.lines_list() for i in [1,-1]])
			Rays = [Ray(i,rays[i][0],rays[i][1],self.num_boundary_comps(rays[i][0]), self.euler_char(rays[i][0]), self.boundary_slopes(rays[i][0]), True, self.is_admissible(rays[i][0]),True, self.simplicial_class(rays[i][0])) for i in range(len(rays))]
			vertices = tuple([(pts_dict[tuple(v)],vector(v)) for v in polyhedron.vertices_list()])
			projected_verts = {}
			Vertices = []
			for i,v in vertices:
				v = vector(v)
				v_proj, coeffs = orthogonal_proj(v, polyhedron.lines_list())
				v_proj = tuple(v_proj)
				boundary_slopes = []
				num_boundary_comps = 0
				euler_char = coeffs[0]*self.euler_char(pts_dict[tuple(v)])
				for j in range(self.manifold.num_cusps()):
					slope = vector((0,0))
					slope += vector(self.boundary_slopes(i)[j])*coeffs[0]
					for k in range(len(polyhedron.lines_list())):
						slope += vector(self.boundary_slopes(rays_dict[tuple(polyhedron.lines_list()[k])])[j])*coeffs[k+1]
					boundary_slopes.append(slope)
					num_boundary_comps += gcd(slope[0],slope[1])
				if v_proj not in projected_verts:
					projected_verts[v_proj] = (num_boundary_comps, euler_char, boundary_slopes)
			polyhedron = Polyhedron(vertices=Matrix(projected_verts.keys()), rays=Matrix(rays_dict.keys()), base_ring=QQ, backend='cdd')
			for i in range(len(polyhedron.vertices_list())):
				v = tuple(polyhedron.vertices_list()[i])
				if v in pts_dict:
					Vertices.append(AdornedVertex(i,pts_dict[v],vector(v),self.num_boundary_comps(pts_dict[v]), self.euler_char(pts_dict[v]), self.boundary_slopes(pts_dict[v]), False, self.is_admissible(pts_dict[v]), True, W.simplicial_class(pts_dict[v])))
				else:
					Vertices.append(AdornedVertex(i, None, vector(v), projected_verts[v][0], projected_verts[v][1], projected_verts[v][2],False, None, True,None))



		else:
			polyhedron = Polyhedron(vertices=Matrix(pts_dict.keys()), base_ring=QQ)
			vertices = tuple([(pts_dict[tuple(v)],vector(v)) for v in polyhedron.vertices_list()])
			if not self.manifold_is_closed:
				Vertices = [AdornedVertex(i,vertices[i][0],vertices[i][1],self.num_boundary_comps(vertices[i][0]), self.euler_char(vertices[i][0]), self.boundary_slopes(vertices[i][0]), False, self.is_admissible(vertices[i][0]), True, self.simplicial_class(vertices[i][0])) for i in range(len(vertices))]
			else:
				Vertices = [AdornedVertex(i,vertices[i][0],vertices[i][1],0, self.euler_char(vertices[i][0]), [], False, self.is_admissible(vertices[i][0]), True, self.simplicial_class(vertices[i][0])) for i in range(len(vertices))]

			Rays = []
		ball = TNormBall(Vertices, Rays, polyhedron)
		if self.num_cusps == 1 and self.betti_number == 1:
			M = self.manifold.copy()
			(p,q) = self.manifold.homological_longitude()
			elem_divs = [div for div in M.homology().elementary_divisors() if div != 0]
			M.dehn_fill((p,q))
			filled_elem_divs = [div for div in M.homology().elementary_divisors() if div != 0]
			b = 1
			for div in elem_divs:
				b *= div
			for div in filled_elem_divs:
				b /= div
			A = self.manifold.alexander_polynomial()
			A_norm = QQ(A.degree() - 1)
			try:
				v=ball.vertices[0]

			except IndexError:
				pass
			if A_norm <= 0 or not A.is_monic():
				self.is_fibered = False
				ball.confirmed = True

			elif len(ball.vertices)==0:
				self.is_fibered = True
				polyhedron = Polyhedron(vertices=[[A_norm],[-A_norm]], base_ring=QQ)
				Vertices = [AdornedVertex(0, None, (-1/A_norm,), b, -A_norm, [(0,1)], False, True, False,None),AdornedVertex(1, None, (1/A_norm,), b, -A_norm, [(0,-1)], False, True, False, None)]
				ball = TNormBall(Vertices, Rays, polyhedron)
				ball.confirmed = True
			elif A_norm == abs(v.euler_char):
				assert v.num_boundary_comps == b
				self.is_fibered = 'unknown'
				ball.confirmed = True
			elif A_norm < abs(v.euler_char):
				M = self.manifold.copy()
				M.dehn_fill(M.homological_longitude())
				Mf = M.filled_triangulation()
				T = regina.Triangulation3(Mf._to_string())
				boo = T.intelligentSimplify()
				ns = regina.NormalSurfaces.enumerate(T,regina.NS_QUAD,regina.NS_VERTEX,regina.NS_ALG_DEFAULT)
				non_trivial = [i for i in range(ns.size()) if ns.surface(i).isOrientable()]
				if len(non_trivial)>1:
					non_trivial = [i for i in non_trivial if ns.surface(i).cutAlong().isConnected()]
				if len(non_trivial) == 0: ## something is wrong!
					print('Warning: failed to confirm that norm ball is correct (error: len(non_trivial)==0, culprit:M={}.'.format(self.manifold.name()))
					ball.confirmed = False
				elif len(non_trivial) >= 1:
					genus = min([(2-regina_to_sage_int(ns.surface(i).eulerChar()))/2 for i in non_trivial])
					if 2*genus - 2 + b == abs(v.euler_char):
						self.is_fibered = False
						ball.confirmed = True
					elif 2*genus - 2 + b == A_norm:
						multiplier = 1/QQ(abs(v.euler_char)/A_norm)
						if self.uses_simplicial_homology == True:
							simplicial_class = multiplier*self.simplicial_class(v.surface_index)
						else:
							simplicial_class = None
						self.is_fibered == True
						polyhedron = Polyhedron(vertices=[[A_norm],[-A_norm]], base_ring=QQ)
						Vertices = [AdornedVertex(0, None, (-1/A_norm,), b, -A_norm, [(0,1)], False, True, False, simplicial_class),AdornedVertex(1, None, (1/A_norm,), b, -A_norm, [(0,-1)], False, True, False, simplicial_class)]
						ball = TNormBall(Vertices, Rays, polyhedron)
						ball.confirmed = True
					else:
						print('Warning: failed to confirm that norm ball is correct (error: 2g-1!=A_norm or abs(X(S)), culprit:M={}.'.format(self.manifold.name()))
						ball.confirmed = False
			else:
				# if we are here, then something is wrong, and the following will cause an error to be thrown.
				assert A_norm <= abs(v.euler_char) # (this should never happen)
		if not self._QUIET:
			print('Done.')
			try:
				sys.stdout.flush()
			except AttributeError:
				pass
		return ball

	@cached_property
	def dual_norm_ball(self):
		B = self.norm_ball
		P = B.polyhedron
		poly_vertices = []
		vertices = []
		for i in range(len(P.faces(P.dim()-1))):
			face = P.faces(P.dim()-1)[i]
			e = face.as_polyhedron().equations_list()[0]
			v = -vector(e[1:])/(e[0])
			poly_vertices.append(v)
			vertices.append(DualVertex(i, v, B.facets(P.dim()-1)[i]))
		P_dual = Polyhedron(vertices=poly_vertices, base_ring=QQ)
		dual_ball = DualNormBall(vertices, P_dual)
		return dual_ball

	@cached_property
	def qtons_info(self):
		if not self._QUIET:
			print('Analyzing quad transversely oriented normal surfaces... ', end='')
			try:
				sys.stdout.flush()
			except AttributeError:
				pass
		qtons_info_dict = dict([(i,{}) for i in range(self.qtons().size())])
		for i in range(self.qtons().size()):
			qtons_info_dict[i]['image_in_H2'] = self.map_to_H2(i)
			qtons_info_dict[i]['euler_char'] = self.euler_char(i)
			qtons_info_dict[i]['boundary_slopes'] = self.boundary_slopes(i)
			qtons_info_dict[i]['num_boundary_comps'] = self.num_boundary_comps(i)
			qtons_info_dict[i]['genus'] = (2-self.euler_char(i)-self.num_boundary_comps(i))/2
			qtons_info_dict[i]['is_norm_minimizing'] = self.is_norm_minimizing(i)
			qtons_info_dict[i]['over_facet'] = self.over_facet(i,as_string=True)
			qtons_info_dict[i]['spinning_slopes'] = self.regina_bdy_slopes(i)
		if not self._QUIET:
			print('Done.')
			try:
				sys.stdout.flush()
			except AttributeError:
				pass
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
			

def over_facet_(v, P):
	if v.is_zero():
		return None
	else:
		ray = Polyhedron(rays=[v]) # ray generated by v
		f_poly = check_subfaces(ray, P) # find the facet that v lies over, as a polyhedron

		# f_poly is the facet we want, but it is in the form of a polyhedron, not a facet of P. So we find the facet of P with the same vertices, and return that.
		for i in range(P.dim()):
			for f in P.faces(i):
				if sorted([vert.vector() for vert in f.vertices()]) == sorted([vert.vector() for vert in f_poly.vertices()]):
					return f

def lies_over_facet(v, P, facet):
	# find the lowest dimensional facet that v lies over
	v_facet = over_facet_(v, P)
	# check if v_facet found above is a subfacet of facet
	if set(v_facet.vertices()).issubset(set(facet.vertices())) and set(v_facet.lines()).issubset(set(facet.lines())):
		return True
	return False

def basis_over_facet(basis, P, facet):  # check if all vectors in basis lie over facet.
	for vec in basis:
		if not lies_over_facet(vec, P, facet):
			return False
	return True

def orthogonal_proj(v, rays):
	coeffs = [QQ(1)]
	for w in rays:
		w = vector(w)
		coeff = QQ(-(v.dot_product(w)/w.dot_product(w)))
		v = v + coeff*w
		coeffs.append(coeff)
	l = lcm([coeff.denominator() for coeff in coeffs])
	if l != 1:
		coeffs = [coeff*l for coeff in coeffs]
		v = v*l
	return v, coeffs

def abs_slopes(slopes):
	return [tuple([abs(int(c)) for c in s]) for s in slopes]

def int_slopes(slopes):
	return [tuple([int(c) for c in s]) for s in slopes]













