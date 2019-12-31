# TNorm

TNorm is a package for computing the Thurston norm unit ball of finite volume orientable hyperbolic 3-manifolds. Currently, tnorm must be installed in Sage, and Sage must have Regina and SnapPy installed. To instal TNorm:

$ sage -pip install tnorm

To run the tnorm graphical user interface app:

$ sage -python -m tnorm.app

To get started:

sage: W=tnorm.load('m130')
sage: B=W.norm_ball
sage: B.vertices
[Vertex 0: represented by (1/2)* S_1,2 at (-1), mapped from surface with index 10,
Vertex 1: represented by (1/2)* S_1,2 at (1), mapped from surface with index 0]
sage: 

In a future release, we plan to remove the dependence on Sage.

Support for hyperbolic 3-manifolds that are not multi-component links in rational homology 3-spheres has been added very recently, and has not been thoroughly tested yet. If you get any results that don't make sense, please email me at william.worden@rice.edu.

TO DO:

* add feature: determine fiberedness of a hyp 3-mfld (and hence knot genus of a knot)
* remove Sage dependence.
* better documentation throughout.
* some optimization for speed is probably still possible.
