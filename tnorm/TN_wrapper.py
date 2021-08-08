#-*-python-*-
from __future__ import print_function
from math import gcd as GCD

import sys

try:
    sys.path.remove('/Applications/Regina.app/Contents/MacOS/python')
except ValueError:
    pass
import regina

import snappy

from tnorm.kernel.simplicial import get_face_map_to_C2, get_quad_map_to_C2, H2_as_subspace_of_C2, qtons_image_in_C2
from tnorm.kernel.boundary import inward_oriented_bdy, outward_oriented_bdy
from tnorm.kernel.euler import solve_lin_gluing_eq, euler_char_
from tnorm.kernel.homology import _map_to_H2, _map_to_H1bdy
from tnorm.kernel.matrices import peripheral_curve_mats, oriented_quads_mat
from tnorm.utilities.regina_helpers import regina_to_sage_int
from tnorm.kernel.embedded import is_embedded, ends_embedded
from tnorm.norm_ball import TNormBall, DualNormBall, NBVertex, Ray, DualVertex
from tnorm.utilities.cached_prop import cached_property
from tnorm.kernel.peripheral import periph_basis_intersections, periph_basis_connected
from tnorm.utilities.sage_types import Matrix, vector, VectorSpace, Polyhedron, RR, QQ
#preparser(False)


class TN_wrapper(object):
	"""
	A class for working with Thurston Norm on H2(M,bdy(M)).

	sage: W = TN_wrapper(SnapPeaTri)

	SnapPeaTri can be any of the following:
	(1) A snappy Manifold or Triangulation
	(2) a regina.SnapPeaTriangulation triangulation
	(3) Any string argument that snappy.Manifold() can take, e.g.,

	sage: W = TN_wrapper('L6a5')
	sage: W = TN_wrapper('/path/to/snappy/file.tri')
	sage: W = TN_wrapper(string)       # where string is the triangulation encoded as a string (i.e., contents of a tri file)
	
	"""
	def __init__(self, manifold, quiet=False, tracker=False, allows_non_admissible=False, force_simplicial_homology=False):
		
		self._triangulation = None
		if isinstance(manifold,str):
			self._manifold = snappy.Manifold(manifold)
		elif isinstance(manifold, regina.engine.SnapPeaTriangulation):
			self._manifold = snappy.Manifold(manifold.snapPea())
			self._triangulation = manifold
		elif isinstance(manifold, regina.engine.Triangulation3):
			self._manifold = snappy.Manifold(manifold.snapPea())
			self._triangulation = manifold
		elif isinstance(manifold, snappy.Manifold):
			self._manifold = manifold
		elif isinstance(manifold, snappy.Triangulation):
			self._manifold = manifold

		for c in self._manifold.cusp_info():
			if c.is_complete == False:
				self._manifold = self._manifold.filled_triangulation()
				break
		
		self._QUIET = quiet
		self._force_simplicial_homology = force_simplicial_homology
		self._num_cusps = self._manifold.num_cusps()
		if self._num_cusps != 0:
			try: 
				L = self._manifold.link()
				self._knows_link_complement = True
			except ValueError:
				self._knows_link_complement = False
			#if self._knows_link_complement and bdy_H1_basis == 'natural':
			#	self._bdy_H1_basis = 'natural'
			#else:
			#	if self._manifold.verify_hyperbolicity()[0]:
			#		self._manifold.set_peripheral_curves('shortest')
			#		self._bdy_H1_basis = 'shortest'

			self._triangulation = regina.SnapPeaTriangulation(self._manifold._to_string())
			self._angle_structure = solve_lin_gluing_eq(self._triangulation)
			self._peripheral_curve_mats = peripheral_curve_mats(self._manifold, self._triangulation)
			self._manifold_is_closed = False


			# check to make sure peripheral basis curves actually give a basis
			check, message = periph_basis_intersections(self)
			assert check, message
			
			# check to make sure each peripheral basis curve is connected (extra trivial loops
			# will mess up Euler char calculations).
			check, message = periph_basis_connected(self)
			assert check, message

		else:
			if self._triangulation == None:
				self._triangulation = regina.Triangulation3(self._manifold._to_string())
				for _ in range(10):
					self._triangulation.intelligentSimplify()
			self._manifold_is_closed = True
			#self._bdy_H1_basis = None
		if not self._triangulation.isOriented():
			self._triangulation.orient()
		self._qtons = False
		self._tkr = False
		self._allows_non_admissible = allows_non_admissible
		self._is_fibered = 'unknown'
		self._betti_number = self._triangulation.homologyH1().rank()

		# qtons memoization caches-----
		self._euler_char = {}
		self._map_to_ball = {}
		self._map_to_H2 = {}
		self._num_boundary_comps = {}
		self._over_facet = {}
		self._is_norm_minimizing = {}
		self._is_admissible = {}
		self._qtons_image_in_C2 = {}
		self._num_H1bdy_comps = {}
		self._has_mixed_bdy = {}
		self._is_embedded = {}
		self._ends_embedded = {}
		self._oriented_quads_mat = {}

		if self._manifold.num_cusps() > 0:
			self._boundary_slopes = {}
			self._spinning_slopes = {}
			self._map_to_H1bdy = {}

		if not self._QUIET:		
			print('Enumerating quad transversely oriented normal surfaces (qtons)... ', end='')
			try:
				sys.stdout.flush()
			except AttributeError:
				pass

		if tracker:
			self._tkr = regina.ProgressTracker()

		# compute transversely oriented normal surfaces
		self._compute_qtons()
		self._qtons = True

		# name the surfaces by their index in the NormalSurfaces list
		for i in range(self.qtons().size()):
			self.qtons().surface(i).setName(str(i))

		if not self._QUIET:
			print('Done.')
			try:
				sys.stdout.flush()
			except AttributeError:
				pass

		if self._betti_number > self._num_cusps:
			self._has_internal_homology = True
		else:
			self._has_internal_homology = False

		# if the manifold has internal homology or force_simplicial_homology==True then we need to use simplicial homology.
		if self._has_internal_homology or self._force_simplicial_homology:
			
			self._uses_simplicial_homology = True

			if not self._QUIET:
				print('computing simplicial homology...',end='')
				try:
					sys.stdout.flush()
				except AttributeError:
					pass

			self._face_map_to_C2 = get_face_map_to_C2(self._triangulation)
			self._quad_map_to_C2 = get_quad_map_to_C2(self._triangulation, self._face_map_to_C2)
			H2_basis_in_C2, P, qtons_image = H2_as_subspace_of_C2(self, self._face_map_to_C2, self._quad_map_to_C2)
			self._project_to_im_del3 = P
			self._qtons_image_in_C2 = {i:qtons_image[i] for i in range(len(qtons_image))}
			assert len(H2_basis_in_C2) == self._betti_number, self._manifold.name()+', force_simplicial_homology={}'.format(self._force_simplicial_homology)
				# raise an error if the qtons do not generate H2. This should only happen for a knot in a 
				# rational homology sphere, with force_simplicial_homology=True.
			I = Matrix.identity(self._betti_number)
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
			self._uses_simplicial_homology = False

	# need to make this work for qtons without index. I.e., a linear combo of vertex normal surfaces.
	

	def manifold(self):
		return self._manifold

	def triangulation(self):
		return self._triangulation

	def num_cusps(self):
		return self._num_cusps

	def knows_link_complement(self):
		return self._knows_link_complement

	#def bdy_H1_basis(self):
	#	return self._bdy_H1_basis

	def allows_non_admissible(self):
		return self._allows_non_admissible

	def betti_number(self):
		return self._betti_number

	def is_fibered(self):
		return self._is_fibered

	def has_internal_homology(self):
		return self._has_internal_homology

	def uses_simplicial_homology(self):
		return self._uses_simplicial_homology

	def manifold_is_closed(self):
		return self._manifold_is_closed


	def _simplicial_map_to_H2(self,qtons):
		if self.uses_simplicial_homology():
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
			if not ind in self._map_to_H2:
				s_H2 = self._map_H2_to_standard_basis*(self._qtons_image_in_C2[ind]-self._project_to_im_del3*self._qtons_image_in_C2[ind])
				self._map_to_H2[ind] = s_H2
			
			return self._map_to_H2[ind]
		else:
			return None

	def _compute_qtons(self):
		if self._tkr:
			if self.allows_non_admissible():
				_ = regina.NormalSurfaces.enumerate(self.triangulation(),regina.NS_ORIENTED_QUAD,regina.NS_IMMERSED_SINGULAR,regina.NS_ALG_DEFAULT,self._tkr)
			else:
				_ = regina.NormalSurfaces.enumerate(self.triangulation(),regina.NS_ORIENTED_QUAD,regina.NS_VERTEX,regina.NS_ALG_DEFAULT,self._tkr)
		else:
			if self.allows_non_admissible():
				_ = regina.NormalSurfaces.enumerate(self.triangulation(),regina.NS_ORIENTED_QUAD,regina.NS_IMMERSED_SINGULAR,regina.NS_ALG_DEFAULT)
			else:
				_ = regina.NormalSurfaces.enumerate(self.triangulation(),regina.NS_ORIENTED_QUAD,regina.NS_VERTEX,regina.NS_ALG_DEFAULT)
		_ = None

	def qtons(self):
		"""
		Enumerate vertex quad transversely oriented normal surfaces for the ideal triangulation self.triangulation()
		This is cached in the Regina triangulation, so calls to this function after the initial call (when TN_wrapper() is 
		instantiated) will be fast. 

		Returns a regina.NormalSurfaces object.
		"""
		if self._qtons == False:
			self._compute_qtons()
			self._qtons == True

		return self._triangulation.lastChild()


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
		-17

		"""
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())
			
		if not ind in self._euler_char:
			s = self.qtons().surface(ind)
			if not self.manifold_is_closed():
				ec = euler_char_(s,self._angle_structure)
			else:
				ec = regina_to_sage_int(s.eulerChar())
			self._euler_char[ind] = ec
		
		return self._euler_char[ind]

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

		if not ind in self._num_boundary_comps:
			if self.manifold_is_closed():
				self._num_boundary_comps[ind] = 0
			else:
				nbc = 0
				pos, neg = self.boundary_slopes(ind).values()
				for slope in pos:
					nbc += gcd(slope[0],slope[1])
				for slope in neg:
					nbc += gcd(slope[0],slope[1])
				self._num_boundary_comps[ind] = nbc

		return self._num_boundary_comps[ind]	

	def boundary_slopes(self,qtons):
		"""
		Return the positive and negative boundary slopes of the given surface in H1(bdy(M)), with respect to the basis of peripheral curves given by SnapPy. 
		Argument can be a surface (regina quad transverely oriented normal surface) or the index of a surface. Returns as a tuple of
		the form (pos_slopes, neg_slopes), where pos_slopes is a list such that pos_slopes[i] is positive boundary slope on cusp i,
		and similarly for neg_slopes.

		sage: W.boundary_slopes(3)    # return positive and negative boudnary of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		([(0, 0), (-3, 1)], [(3, -1), (0, 0)])
		"""
		if self.manifold_is_closed():
			return []
		else:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if not ind in self._boundary_slopes:
				s = self.qtons().surface(ind)
				out_bdy = outward_oriented_bdy(s, self._peripheral_curve_mats, self.oriented_quads_mat(s), False)
				in_bdy = inward_oriented_bdy(s, self._peripheral_curve_mats, self.oriented_quads_mat(s), False)
				self._boundary_slopes[ind] = {'outward':out_bdy, 'inward':in_bdy}
			
			return self._boundary_slopes[ind]

	def oriented_quads_mat(self, qtons):
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())

		if not ind in self._oriented_quads_mat:
			s = self.qtons().surface(ind)
			self._oriented_quads_mat[ind] = oriented_quads_mat(s)
		
		return self._oriented_quads_mat[ind]

	def is_embedded(self,qtons):
		"""Return True if the given qtons can be normally isotoped to be embedded.
		"""
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())

		if not ind in self._is_embedded:
			s = self.qtons().surface(ind)
			self._is_embedded[ind] = is_embedded(s,self)
		
		return self._is_embedded[ind]

	def ends_embedded(self,qtons):
		"""Return True if the ends of the given qtons can be normally isotoped to be embedded.
		"""
		if self.manifold_is_closed():
			return True

		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())

		if not ind in self._ends_embedded:
			s = self.qtons().surface(ind)
			self._ends_embedded[ind] = ends_embedded(s,self)
		
		return self._ends_embedded[ind]		

	def has_mixed_bdy(self, qtons):
		if self.manifold_is_closed():
			return False
		else:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if not ind in self._has_mixed_bdy:
				pos_bdy, neg_bdy = self.boundary_slopes(ind).values()
				if Matrix(pos_bdy).norm() == 0 or Matrix(neg_bdy).norm() == 0:
					self._has_mixed_bdy[ind] = False
				else:
					self._has_mixed_bdy[ind] = True

			return self._has_mixed_bdy[ind]


	def spinning_slopes(self,qtons):
		"""
		Return the positive and negative boundary slopes of the given surface in H1(bdy(M)), with respect to the basis of peripheral curves given by SnapPy.
		Argument can be a surface (regina quad transverely oriented normal surface) or the index of a surface.

		sage: W._map_to_H1bdy(3)    # return image in H1(bdy(M)) of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		[(-1, 0), (0, 1), (-1, 0)]
		"""
		if self.manifold_is_closed():
			return []
		else:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if not ind in self._spinning_slopes:
				s = self.qtons().surface(ind)
				out_bdy = outward_oriented_bdy(s, self._peripheral_curve_mats, self.oriented_quads_mat(s), True)
				in_bdy = inward_oriented_bdy(s, self._peripheral_curve_mats, self.oriented_quads_mat(s), True)
				self._spinning_slopes[ind] = out_bdy, in_bdy
			
			return self._spinning_slopes[ind]



	def H1bdy_slopes(self, qtons):
		
		if self.manifold_is_closed():
			return []
		else:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if not ind in self._map_to_H1bdy:
				s = self.qtons().surface(ind)
				slopes = _map_to_H1bdy(s, self)
				self._map_to_H1bdy[ind] = slopes
			
			return self._map_to_H1bdy[ind]	

	def num_H1bdy_comps(self, qtons):
		try:
			ind = int(qtons)
			s = self.qtons().surface(ind)
		except TypeError:
			ind = int(qtons.name())
			s = qtons

		if not ind in self._num_H1bdy_comps:
			nbc = 0
			slopes = self.H1bdy_slopes(ind)
			for slope in slopes:
				nbc += gcd(slope[0],slope[1])
			self._num_H1bdy_comps[ind] = nbc

		return self._num_H1bdy_comps[ind]	


	def simplicial_class(self,qtons):
		""" Return the vector in the dimension 2 chain group C2 corresponding to the qtons surface. The i^th coordinate
			of this vector corresponds to the number of triangles of index i, as indexed by Regina, with orientation positive
			if the transverse orientation of the surface points toward the tetrahedron triangle(i).front(). 
		"""
		if not self.uses_simplicial_homology():
			return None
		else:
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if not ind in self._qtons_image_in_C2:
				c = qtons_image_in_C2(self,self._quad_map_to_C2)
				self._qtons_image_in_C2[ind] = c

			return self._qtons_image_in_C2[ind]

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
			f = over_facet_(v,self.norm_ball.polyhedron())
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

		if not ind in self._is_norm_minimizing:
			v = self.map_to_H2(ind)
			if not self.norm_ball.polyhedron().interior_contains(v):
				self._is_norm_minimizing[ind] = True
				return True
			else:
				self._is_norm_minimizing[ind] = False

		return self._is_norm_minimizing[ind]

	def is_admissible(self, qtons):
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())

		if not ind in self._is_admissible:
			s = self.qtons().surface(ind)
			mat = oriented_quads_mat(s)
			for row in mat:
				nonzero = [i for i in range(0,6,2) if ( row[i]!=0 or row[i+1]!=0 )]
				if len(nonzero) > 1:
					self._is_admissible[ind] = False
					return False
			self._is_admissible[ind] = True

		return self._is_admissible[ind]

	def map_to_ball(self, qtons):
		"""
		Return the image of the given surface in H2(M, bdy(M)), then divide by Euler characteristic. Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		"""
		try:
			ind = int(qtons)
		except TypeError:
			ind = int(qtons.name())

		if not ind in self._map_to_ball:
			image_in_H2 = self.map_to_H2(ind)
			if self.euler_char(ind) < 0:
				mtb = tuple([QQ(i)/(-self.euler_char(ind)) for i in image_in_H2])
				self._map_to_ball[ind] = mtb
				return mtb
			else:
				self._map_to_ball[ind] = float('inf')

		return self._map_to_ball[ind]


	def map_to_H2(self, qtons):
		"""
		Return the image of the given surface in H2(M, bdy(M)). Argument can be a surface (regina quad transversely oriented
		normal surface) or the index of a surface.

		sage: W.mapToH2(3)    # return the image in H2(M,bdy(M)) of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		(0,1,0)

		"""
		if not self.uses_simplicial_homology():
			try:
				ind = int(qtons)
			except TypeError:
				ind = int(qtons.name())
	
			if not ind in self._map_to_H2:
				s = self.qtons().surface(ind)
				self._map_to_H2[ind] = _map_to_H2(s, self)
			
			return self._map_to_H2[ind]
		else:
			return self._simplicial_map_to_H2(qtons)
		
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
				good_ends = [rep for rep in reps if self.ends_embedded(rep)]
				if len(good_ends) != 0:
					reps = good_ends # if possible, choose from reps with embedded ends
				embedded = [rep for rep in reps if self.is_embedded(rep)]
				if len(embedded) != 0:
					reps = embedded # if possible, choose from reps that are embedded
				reps = sorted(reps,key=lambda x: self.num_boundary_comps(x))
				points[nzd] = reps[0] # get the index of the rep with the minimal number of bdy components

		return points, rays

	@cached_property
	def norm_ball(self):
		"""
		Return the Thurston norm ball.
		"""
		if not self._QUIET:
			print('Computing Thurston norm unit ball... ', end='')
			try:
				sys.stdout.flush()
			except AttributeError:
				pass
		pts_dict, rays_dict = self._norm_ball_points

		### If M is not hyperbolic then rays_dict should be non-empty, and the norm ball should be non-compact. 
		### The below attempts to compute the non-compact norm ball, but it may be wrong. This is because 
		### (it seems) some surfaces may not be realized as spun normal surfaces.

		if len(rays_dict) != 0:
			V=VectorSpace(RR,self.betti_number())
			ray_span=V.span([V(r) for r in rays_dict])
			keys = [p for p in pts_dict]
			for p in keys:
				if V(p) in ray_span:
					_ = pts_dict.pop(p)


			polyhedron = Polyhedron(vertices=Matrix(pts_dict.keys()), rays=Matrix(rays_dict.keys()), base_ring=QQ, backend='cdd')
			rays = tuple([(rays_dict[tuple(i*vector(v))],i*vector(v)) for v in polyhedron.lines_list() for i in [1,-1]])
			Rays = [Ray(i, rays[i][0], rays[i][1], self) for i in range(len(rays))]
			vertices = tuple([(pts_dict[tuple(v)],vector(v)) for v in polyhedron.vertices_list() if not vector(v).is_zero()])
			Vertices = [NBVertex(i, vertices[i][0], vertices[i][1], self) for i in range(len(vertices))]

		else:
			polyhedron = Polyhedron(vertices=Matrix(pts_dict.keys()), base_ring=QQ)
			vertices = tuple([(pts_dict[tuple(v)],vector(v)) for v in polyhedron.vertices_list()])
			
			Vertices = [NBVertex(i,vertices[i][0],vertices[i][1], self) for i in range(len(vertices))]
			Rays = []

		ball = TNormBall(Vertices, Rays, polyhedron)
		if self.num_cusps() == 1 and self.betti_number() == 1:
			M = self.manifold().copy()
			(p,q) = self.manifold().homological_longitude()
			elem_divs = [div for div in M.homology().elementary_divisors() if div != 0]
			M.dehn_fill((p,q))
			filled_elem_divs = [div for div in M.homology().elementary_divisors() if div != 0]
			b = 1
			for div in elem_divs:
				b *= div
			for div in filled_elem_divs:
				b /= div
			A = self.manifold().alexander_polynomial()
			A_norm = QQ(A.degree() - 1)
			try:
				v=ball.vertices()[0]

			except IndexError:
				pass

			### below needs to be fixed. Currently does not account for virtual fibers.
			if A_norm <= 0 or not A.is_monic():
				self._is_fibered = False
				ball._confirmed = True

			elif len(ball.vertices())==0:
				self._is_fibered = True
				polyhedron = Polyhedron(vertices=[[A_norm],[-A_norm]], base_ring=QQ)
				Vertices = [NBVertex(0, None, (-1/A_norm,), self, (b, -A_norm, {'outward':[(0,1)],'inward':[(0,0)]})), NBVertex(1, None, (1/A_norm,), self, (b, -A_norm, {'outward':[(0,-1)],'inward':[(0,0)]}))]
				ball = TNormBall(Vertices, Rays, polyhedron)
				ball._confirmed = True
			elif A_norm == abs(v.euler_char()):
				assert self.num_H1bdy_comps(v.qtons_index()) == b
				self._is_fibered = 'unknown'
				ball._confirmed = True
			elif A_norm < abs(v.euler_char()):
				#M = self.manifold().copy()
				#M.dehn_fill(M.homological_longitude())
				Mf = M.filled_triangulation()
				T = regina.Triangulation3(Mf._to_string())
				boo = T.intelligentSimplify()
				ns = regina.NormalSurfaces.enumerate(T,regina.NS_QUAD,regina.NS_VERTEX,regina.NS_ALG_DEFAULT)
				non_trivial = [i for i in range(ns.size()) if ns.surface(i).isOrientable()]
				if len(non_trivial)>1:
					non_trivial = [i for i in non_trivial if ns.surface(i).cutAlong().isConnected()]
				if len(non_trivial) == 0: ## something is wrong!
					print('Warning: failed to confirm that norm ball is correct (error: len(non_trivial)==0, culprit:M={}.'.format(self.manifold().name()))
					ball._confirmed = False
				elif len(non_trivial) >= 1:
					genus = min([(2-regina_to_sage_int(ns.surface(i).eulerChar()))/2 for i in non_trivial])
					if 2*genus - 2 + b == abs(v.euler_char()):
						self._is_fibered = False
						ball._confirmed = True
					elif 2*genus - 2 + b == A_norm:
						multiplier = 1/QQ(abs(v.euler_char())/A_norm)
						if self.uses_simplicial_homology() == True:
							simplicial_class = multiplier*self.simplicial_class(v.qtons_index())
						else:
							simplicial_class = None
						self._is_fibered == True
						polyhedron = Polyhedron(vertices=[[A_norm],[-A_norm]], base_ring=QQ)
						Vertices = [NBVertex(0, None, (-1/A_norm,), self, (b, -A_norm, {'outward':[(0,1)],'inward':[(0,0)]})),NBVertex(1, None, (1/A_norm,), self, (b, -A_norm, {'outward':[(0,-1)],'inward':[(0,0)]}))]
						ball = TNormBall(Vertices, Rays, polyhedron)
						ball._confirmed = True
					else:
						print('Warning: failed to confirm that norm ball is correct (error: 2g-1!=A_norm or abs(X(S)), culprit:M={}.'.format(self.manifold().name()))
						ball._confirmed = False
			else:
				# if we are here, then something is wrong, and the following will cause an error to be thrown.
				assert A_norm <= abs(v.euler_char()) # (this should never happen)
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
		P = B.polyhedron()
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
			qtons_info_dict[i]['is_embedded'] = self.is_embedded(i)
			qtons_info_dict[i]['num_boundary_comps'] = self.num_boundary_comps(i)
			qtons_info_dict[i]['genus'] = self.genus(i)
			qtons_info_dict[i]['is_norm_minimizing'] = self.is_norm_minimizing(i)
			qtons_info_dict[i]['over_facet'] = self.over_facet(i,as_string=True)
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


def lcm(list_of_ints):
    L = list_of_ints
    LCM = 1
    for k in L:
        LCM = LCM*k//gcd(LCM,k)
    return LCM

def gcd(a,b):
    return abs(GCD(a,b))








