Lset-opt is a level-set package customized for topological optimization.

Lset-opt tracks the contours of shapes using a custom 2D level-set method. The level-set function, phi, defines the contours of the shapes on the grid. The contour then defines the filling fraction, p, of each cell on the grid.

The primary functionality of lset-opt is the ability to update the contours based on a desired change in p. This is accomplished by updating phi (the level-set function) to reflect the changes in p which are near the edge of shapes. The ability to generate new shapes (e.g. small islands or holes) is included as well.


Regularization
--------------

Phi, the level-set function, is regularized simply by specifying a level-set function which defines the shape boundary on its zero level-set (where it crosses zero). Critically, regularization is guaranteed to preserve all shape boundaries; certain elements in phi are simply resized in an (amateurish) attempt to increase numerical accuracy. Phi can be regularized with arbitrary frequency, although back-to-back regularizations produce identical phi as that produced by a single regularization.

A single exception occurs when a value of phi equals 0, in this case that value is set to the smallest positive number available. This is done so that every boundary point is always defined by exactly two values of phi.


Boundary definition
-------------------

Initialization is able to exactly preserve shape boundaries because boundary points are simply defined as the zero-crossing between adjacent grid points whose values of phi are of opposite sign. Only horizontally- or vertically-adjacent point pairs are considered, and the zero-crossing is determined by linear interpolation. 


Fill-fraction
-------------

The fill-fraction, p, at every grid point is determined from the boundary points and the polarity of phi. The values of p may range from -1 to +1, inclusive. A value of -1 indicates that the cell is completely filled with material A, while +1 indicates complete filling by material B. 

Lset-opt defines each cell as the box centered at a grid point, with length and height of 2.0 (the grid spacing is 1.0). This causes the cells to overlap, or smear with one another, which seems more useful than non-overlapping cells. Specifically, overlapping cells provide more information as to which direction the boundary should move in.

    The sign of phi determines the material at the center of the cell, which will fill a rectangle within the cell. The remainder of the cell is filled with the other material, and the fill-fraction is given by:

    p = 0.5 * sign(phi) * (width-1 + height-1).
    

The width of the inner rectangle is determined by the boundary points to the left and right of the grid point. If no boundary point exists on a side, then we assume that the inner rectangle extends a distance of 1.0 in that direction. The height of the rectangle is calculated in an equivalent way.


Dynamic shapes
--------------

Lset-opt allows one to update the topology of the grid by specifying a proposed change in p, dp. P is updated indirectly by changing phi, which guarantees that p is still based on a two-material topology. At the same time, this will of course cause considerable difference between the new p and p + dp, which the user should expect.


Island creation
---------------

The definition of p allows for the creation of arbitrarily small "islands" of material. Although, numerous arbitrarily thin islands are possible, lset-opt only allows islands of one-grid point to be created. The creation of possibly multiple islands occurs automatically when updating p, although this feature can also be turned off.
