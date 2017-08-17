lat.Const=[1]                      // lattice constant, 1x1 real
lat.Primitive=[1]                  // Primitive vectors, (3x3/2x2/1x1)
lat.Sublatt=[0;0.1]                    // sublattice position, (nx3/nx2/nx1/)
lat.Order=[2]                     // Nearest Neighbor Order, 1x1 integer

============= PiLib Variable =============
lat.recip_vec, @full, the reciprocal lattice vectors
ORDER=    0, SIZE=[    1,    1], TYPE=REAL

           1

    6.283185

============= PiLib Variable =============
lat.surr_site(1), @full, surrouding sites [order, dist, sublatt, n1, n2, n3, x, y, z]
ORDER=    0, SIZE=[    3,    9], TYPE=REAL

           1           2           3           4           5           6           7           8           9

    0.000000    0.000000    1.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
    1.000000    0.100000    2.000000    0.000000    0.000000    0.000000    0.100000    0.000000    0.000000
    2.000000    0.900000    2.000000   -1.000000    0.000000    0.000000   -0.900000    0.000000    0.000000

============= PiLib Variable =============
lat.surr_site(2), @full, surrouding sites [order, dist, sublatt, n1, n2, n3, x, y, z]
ORDER=    0, SIZE=[    3,    9], TYPE=REAL

           1           2           3           4           5           6           7           8           9

    0.000000    0.000000    2.000000    0.000000    0.000000    0.000000    0.100000    0.000000    0.000000
    1.000000    0.100000    1.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
    2.000000    0.900000    1.000000    1.000000    0.000000    0.000000    1.000000    0.000000    0.000000
