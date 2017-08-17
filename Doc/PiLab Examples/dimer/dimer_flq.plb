flq.Frequency=[10]                  // field frequency, 1x1 real
flq.Order=[3]                      // order of photon process, 1x1, int
flq.Amplitude=[7.2]                  // AC amplitude, 1x1 / 1x2 / 1x3, real
flq.Phase=[-%pi/2]                      // AC phase 1x1 / 1x2 / 1x3, real

============= PiLib Variable =============
flq.state_info, [state_label, order, site, n, l, SubOrb]
ORDER=    0, SIZE=[   28,    6], TYPE=INTEGER

       1       2       3       4       5       6

       1      -3       1       1       0       1
       2      -3       1       1       0       2
       3      -3       2       1       0       1
       4      -3       2       1       0       2
       5      -2       1       1       0       1
       6      -2       1       1       0       2
       7      -2       2       1       0       1
       8      -2       2       1       0       2
       9      -1       1       1       0       1
      10      -1       1       1       0       2
      11      -1       2       1       0       1
      12      -1       2       1       0       2
      13       0       1       1       0       1
      14       0       1       1       0       2
      15       0       2       1       0       1
      16       0       2       1       0       2
      17       1       1       1       0       1
      18       1       1       1       0       2
      19       1       2       1       0       1
      20       1       2       1       0       2
      21       2       1       1       0       1
      22       2       1       1       0       2
      23       2       2       1       0       1
      24       2       2       1       0       2
      25       3       1       1       0       1
      26       3       1       1       0       2
      27       3       2       1       0       1
      28       3       2       1       0       2

============= PiLib Variable =============
flq.H_onsite(1) (t-sp), Floquet H_onsite of order 0
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000

============= PiLib Variable =============
flq.H_onsite(2) (t-sp), Floquet H_onsite of order 1
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000

============= PiLib Variable =============
flq.H_onsite(3) (t-sp), Floquet H_onsite of order 2
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000

============= PiLib Variable =============
flq.H_onsite(4) (t-sp), Floquet H_onsite of order 3
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000

============= PiLib Variable =============
flq.hop_size, size of flq.hop_mat, [order+1,sublatt,hop_mat_size]
ORDER=    0, SIZE=[    8,    5], TYPE=INTEGER

       1       2       3       4       5

       1       1       4       4       2
       1       2       4       4       2
       2       1       4       4       2
       2       2       4       4       2
       3       1       4       4       2
       3       2       4       4       2
       4       1       4       4       2
       4       2       4       4       2

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,1) (t-sp), Floquet hop_mat(1)(:,:,1) of order 0
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.262362  0.000000
       2       4    0.262362  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,2) (t-sp), Floquet hop_mat(1)(:,:,2) of order 0
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.256961  0.000000
       2       4    0.256961  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,1) (t-sp), Floquet hop_mat(2)(:,:,1) of order 0
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.262362  0.000000
       2       4    0.262362  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,2) (t-sp), Floquet hop_mat(2)(:,:,2) of order 0
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.256961  0.000000
       2       4    0.256961  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,1) (t-sp), Floquet hop_mat(1)(:,:,1) of order 1
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.101151  0.000000
       2       4    0.101151  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,2) (t-sp), Floquet hop_mat(1)(:,:,2) of order 1
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.159495  0.000000
       2       4    0.159495  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,1) (t-sp), Floquet hop_mat(2)(:,:,1) of order 1
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3   -0.101151 -0.000000
       2       4   -0.101151 -0.000000

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,2) (t-sp), Floquet hop_mat(2)(:,:,2) of order 1
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3   -0.159495 -0.000000
       2       4   -0.159495 -0.000000

============= PiLib Variable =============
flq.hop_mat(3)(1)(:,:,1) (t-sp), Floquet hop_mat(1)(:,:,1) of order 2
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.018614  0.000000
       2       4    0.018614  0.000000

============= PiLib Variable =============
flq.hop_mat(3)(1)(:,:,2) (t-sp), Floquet hop_mat(1)(:,:,2) of order 2
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3   -0.306188 -0.000000
       2       4   -0.306188 -0.000000

============= PiLib Variable =============
flq.hop_mat(3)(2)(:,:,1) (t-sp), Floquet hop_mat(2)(:,:,1) of order 2
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.018614  0.000000
       2       4    0.018614  0.000000

============= PiLib Variable =============
flq.hop_mat(3)(2)(:,:,2) (t-sp), Floquet hop_mat(2)(:,:,2) of order 2
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3   -0.306188 -0.000000
       2       4   -0.306188 -0.000000

============= PiLib Variable =============
flq.hop_mat(4)(1)(:,:,1) (t-sp), Floquet hop_mat(1)(:,:,1) of order 3
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.002258  0.000000
       2       4    0.002258  0.000000

============= PiLib Variable =============
flq.hop_mat(4)(1)(:,:,2) (t-sp), Floquet hop_mat(1)(:,:,2) of order 3
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3    0.029510  0.000000
       2       4    0.029510  0.000000

============= PiLib Variable =============
flq.hop_mat(4)(2)(:,:,1) (t-sp), Floquet hop_mat(2)(:,:,1) of order 3
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3   -0.002258 -0.000000
       2       4   -0.002258 -0.000000

============= PiLib Variable =============
flq.hop_mat(4)(2)(:,:,2) (t-sp), Floquet hop_mat(2)(:,:,2) of order 3
ORDER=    0, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

       4       4    0.000000  0.000000
       1       3   -0.029510 -0.000000
       2       4   -0.029510 -0.000000
