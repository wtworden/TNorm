#-*-python-*-


import regina
import snappy
import sys
import time


from tnorm.eulerChar import *
from tnorm.boundary import *
from tnorm.HashableDict import *
from tnorm.homology import *
from tnorm.regina_helpers import *
from tnorm.normBall import *
from tnorm.sage_types import *
from tnorm.matrices import *


#preparser(False)


def cachedproperty(func):
    """ Used on methods to convert them to methods that replace themselves\
        with their return value once they are called. """

    def cache(*args):
        self = args[0] # Reference to the class who owns the method
        funcname = func.__name__
        ret_value = func(self)
        setattr(self, funcname, ret_value) # Replace the function with its value
        return ret_value # Return the result of the function

    return property(cache)


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
	def __init__(self,SnapPeaTri,QTONS=None,tracker=False,quiet=False,allowsNonAdmissible=False):
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
		self._normBallPoints = None
		self._QTONS = QTONS
		self._index_bound = None
		self._tkr = False
		self.allowsNonAdmissible = allowsNonAdmissible
		if self._QTONS == None:
			if not quiet:			
				print("enumerating oriented spun normal surfaces. \nThis could take awhile if the triangulation is large!")
			if tracker:
				self._tkr = regina.ProgressTracker()
				if not quiet:
					print('Warning: we are using a progress tracker. \nThe sage prompt will be returned imediately, but computations are\n continuing in the background. Pass W.showProgress() to check if finished.')
			self._Oriented_QTONS_List()



	def _Oriented_QTONS_List(self):
		"""
		Enumerate vertex oriented quad normal surfaces for the ideal triangulation self.triangulation()
		This is cached in self._QTONS, so calls to this function after the initial call (when TN_wrapper() is 
		instantiated) will be fast. 

		Returns a regina.NormalSurfaces object.
		"""
		if self._QTONS == None:
			if self._tkr:
				if self.allowsNonAdmissible:
					self._QTONS = regina.NormalSurfaces.enumerate(self.triangulation,regina.NS_ORIENTED_QUAD,regina.NS_IMMERSED_SINGULAR,regina.NS_ALG_DEFAULT,self._tkr)
				else:
					self._QTONS = regina.NormalSurfaces.enumerate(self.triangulation,regina.NS_ORIENTED_QUAD,regina.NS_VERTEX,regina.NS_ALG_DEFAULT,self._tkr)
			else:
				if self.allowsNonAdmissible:
					self._QTONS = regina.NormalSurfaces.enumerate(self.triangulation,regina.NS_ORIENTED_QUAD,regina.NS_IMMERSED_SINGULAR,regina.NS_ALG_DEFAULT)
				else:
					self._QTONS = regina.NormalSurfaces.enumerate(self.triangulation,regina.NS_ORIENTED_QUAD,regina.NS_VERTEX,regina.NS_ALG_DEFAULT)
		else:
			pass

	def showProgress(self):
		print(self._tkr.percent())

	def orientedNormalSurfaceList(self):
		self._Oriented_QTONS_List()
		return self._QTONS

	def orientedNormalSurface(self,index):
		return self.orientedNormalSurfaceList().surface(int(index))

	def eulerChar(self,sn_surface):
		"""
		Return the Euler characteristic of the spun normal surface. 
		Argument can be a surface (Regina oriented spun normal surface) or the index of a surface.

		sage: W.eulerChar(3)    # return Euler characteristic of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		-1

		"""
		if isinstance(sn_surface,int) or isinstance(sn_surface, Integer):
			return eulerChar(self._QTONS.surface(int(sn_surface)),self._angle_structure)
		else:
			return eulerChar(sn_surface,self._angle_structure)

	def boundarySlopes(self,sn_surface):
		"""
		Return the boundary slopes of the spun normal surface. This will be different than the return value
		of Regina's method boundaryIntersections(), because Regina returns the intersection numbers,
		not the slope (i.e., Regina's method will be the negative inverse of these slopes). 
		Argument can be a surface (Regina oriented spun normal surface) or the index of a surface.

		sage: W.boundarySlopes(3)    # return boundary slopes of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		[(-1, 0), (0, 1), (-1, 0)]

		"""
		if isinstance(sn_surface,int) or isinstance(sn_surface, Integer):
			return boundarySlopes(self._QTONS.surface(int(sn_surface)),self.manifold)
		else:
			return boundarySlopes(sn_surface,self.manifold)

	def _map_to_H1bdy(self,sn_surface):
		"""
		Return the image of the given surface in H1(\partial M). Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		sage: W._map_to_H1bdy(3)    # return image in H1(\partial M) of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		[(-1, 0), (0, 1), (-1, 0)]

		"""
		if isinstance(sn_surface,int) or isinstance(sn_surface, Integer):
			return Q_to_H1bdy(self._QTONS.surface(int(sn_surface)),self.manifold)
		else:
			return Q_to_H1bdy(sn_surface,self.manifold)
	
	def mapToH2(self,sn_surface):
		"""
		Return the image of the given surface in H2(M, \partial M). Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		sage: W.mapToH2(3)    # return the image in H2(M,\partial M) of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		(0,1,0)

		"""
		if isinstance(sn_surface,int) or isinstance(sn_surface, Integer):
			return homology_map(self._QTONS.surface(int(sn_surface)),self.manifold)
		else:
			return homology_map(sn_surface,self.manifold)	

	def mapToBall(self,sn_surface):
		"""
		Return the image of the given surface in H2(M, \partial M), then divide by Euler characteristic. Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		"""
		imageInH2 = self.mapToH2(sn_surface)
		if self.eulerChar(sn_surface) < 0:
			return tuple([i/(-self.eulerChar(sn_surface)) for i in imageInH2])
		else:
			return None


	def genus(self,sn_surface):
		"""
		Return the genus. Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		sage: W.genus(3)    # return genus of the surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		1

		"""
		return (2-self.eulerChar(sn_surface)-self.numBoundaryComps(sn_surface))/2
	
	def numBoundaryComps(self,sn_surface):
		"""
		Return the number of boundary components. Argument can be a surface (regina oriented
		spun normal surface) or the index of a surface.

		sage: W.numBoundaryComps(3)    # return num of boundary components of surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		2

		"""	
		if isinstance(sn_surface,int) or isinstance(sn_surface, Integer):
			return numBoundaryComps(self._QTONS.surface(int(sn_surface)))
		else:
			return numBoundaryComps(sn_surface)

	def numBoundaryCompsH1(self,sn_surface):
		"""
		Return the number of boundary components after mapping to H1(bdy m). This will be the number
		of boundary components after cancelling oppositely oriented pairs (i.e., a 3-puntured sphere
		will become a once-punctured torus if two of the punctures are on the cusps and boundary curves
		with opposite orientations). Argument can be a surface (regina oriented spun normal surface) or 
		the index of a surface.

		sage: W.numBoundaryComps(3)    # return num of boundary components of surface with index 3, i.e, W.OrientedNormalSurfaceList.surface(3).
		2

		"""	
		if isinstance(sn_surface,int) or isinstance(sn_surface, Integer):
			S = self._QTONS.surface(int(sn_surface))
		else:
			S = sn_surface

		numBoundaryComps = 0
		b = self._map_to_H1bdy(S)
		for slope in b:
			numBoundaryComps += gcd(slope[0],slope[1])
		return numBoundaryComps

		
	def _image_in_H2M(self):
		image = []
		for i in range(self._index_bound):
			pt = self.mapToH2(i)
			if not pt.is_zero():
				ec = self.eulerChar(i)
				if (pt,ec) not in [tup[1:] for tup in image]:
					image.append((i,pt,ec))
		return image

#	@cachedproperty
	def normBallPoints(self):
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
		return HashableDict(image)


	def normBall(self,index_bound=None):
		"""
		Return the norm ball as a sage Polyhedron. Some things you can do:

		sage: B=W.normBall()
		sage: B.faces(2)        # list the 2-dimensional faces
		sage: B.faces(0)        # list the 0-dimensional faces (vertices)

		sage: F = B.faces(2)
		sage: F[0].vertices()   #return vertices of the 0 face

		sage: B.plot()          # plot of B. I imagine this only works if B is 1,2, or 3 dimensional.

		The Polyhedron class has a lot of methods, and I have not explored many of them. There is 
		probably a lot more that would be useful. Use tab completion to explore!
		"""
		if index_bound == None:
			self._index_bound = self._QTONS.size()
		else:
			self._index_bound = index_bound
		pts_dict = self.normBallPoints()
		poly = Polyhedron(vertices=Matrix(pts_dict.keys()))
		verts = tuple([(pts_dict[tuple(v)],vector(v)) for v in poly.vertices_list()])
		vertices = VerticesList(tuple([Vertex(i,verts[i][0],verts[i][1],self.numBoundaryCompsH1(verts[i][0]),self.eulerChar(verts[i][0]),self.boundarySlopes(verts[i][0]),self._map_to_H1bdy(verts[i][0])) for i in range(len(verts))]))
		return NormBall(vertices, poly)


























