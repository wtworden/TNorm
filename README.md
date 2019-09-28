# TNorm

TNorm is a package for computing the Thurston norm unit ball of hyperbolic links in homology 3-spheres. Currently, tnorm must be installed in Sage, and Sage must have Regina and SnapPy installed. To instal TNorm:

$ sage -pip install tnorm

To run the tnorm graphical user interface app:

$ sage -python -m tnorm.app

In a future release, we plan to remove the dependence on Sage.

Currently, if M is a fibered knot complement, tnorm might not return the correct norm ball. We plan to fix this issue in a future release. Additionally, if M is not hyperbolic, the Thurston norm is only a semi-norm, and the unit norm ball is an infinite polyhedron. In this case, tnorm will do its best, but may not return the correct norm ball.