Lset-opt is a level-set package customized for topological optimization.

Description
===========

Lset-opt tracks the contours of shapes using a custom 2D level-set method. The level-set function, phi, defines the contours of the shapes on the grid. The contour then defines the filling fraction, p, of each cell on the grid.

The primary functionality of lset-opt is the ability to update the contours based on a desired change in p. This is accomplished by updating phi (the level-set function) to reflect the changes in p which are near the edge of shapes. The ability to generate new shapes (e.g. small islands or holes) is included as well.


Shape initialization
====================

To determine the initial topology a 
