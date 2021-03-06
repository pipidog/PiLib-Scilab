lat.Const=[4.54]                     // lattice constant, 1x1 real
lat.Primitive=...                 // Primitive vectors, (3x3/2x2/1x1)
[1/2,1/2,0;1/2,0,1/2;0,1/2,1/2]
lat.Sublatt=[0,0,0;1/2,1/2,1/2]   // sublattice position, (nx3/nx2/nx1/)
lat.Order=[1]                     // Nearest Neighbor Order, 1x1 integer

============= PiLib Variable =============
lat.recip_vec, @full, the reciprocal lattice vectors
ORDER=    0, SIZE=[    3,    3], TYPE=REAL

           1           2           3

    1.383962    1.383962   -1.383962
    1.383962   -1.383962    1.383962
   -1.383962    1.383962    1.383962

============= PiLib Variable =============
lat.surr_site(1), @full, surrouding sites [order, dist, sublatt, n1, n2, n3, x, y, z]
ORDER=    0, SIZE=[    7,    9], TYPE=REAL

           1           2           3           4           5           6           7           8           9

    0.000000    0.000000    1.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
    1.000000    2.270000    2.000000   -1.000000   -1.000000    0.000000   -2.270000    0.000000    0.000000
    1.000000    2.270000    2.000000   -1.000000    0.000000   -1.000000    0.000000   -2.270000    0.000000
    1.000000    2.270000    2.000000   -1.000000    0.000000    0.000000    0.000000    0.000000    2.270000
    1.000000    2.270000    2.000000    0.000000   -1.000000   -1.000000    0.000000    0.000000   -2.270000
    1.000000    2.270000    2.000000    0.000000   -1.000000    0.000000    0.000000    2.270000    0.000000
    1.000000    2.270000    2.000000    0.000000    0.000000   -1.000000    2.270000    0.000000    0.000000

============= PiLib Variable =============
lat.surr_site(2), @full, surrouding sites [order, dist, sublatt, n1, n2, n3, x, y, z]
ORDER=    0, SIZE=[    7,    9], TYPE=REAL

           1           2           3           4           5           6           7           8           9

    0.000000    0.000000    2.000000    0.000000    0.000000    0.000000    2.270000    2.270000    2.270000
    1.000000    2.270000    1.000000    0.000000    0.000000    1.000000    0.000000    2.270000    2.270000
    1.000000    2.270000    1.000000    0.000000    1.000000    0.000000    2.270000    0.000000    2.270000
    1.000000    2.270000    1.000000    0.000000    1.000000    1.000000    2.270000    2.270000    4.540000
    1.000000    2.270000    1.000000    1.000000    0.000000    0.000000    2.270000    2.270000    0.000000
    1.000000    2.270000    1.000000    1.000000    0.000000    1.000000    2.270000    4.540000    2.270000
    1.000000    2.270000    1.000000    1.000000    1.000000    0.000000    4.540000    2.270000    2.270000
