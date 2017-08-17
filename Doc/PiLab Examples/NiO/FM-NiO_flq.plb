flq.Frequency=[15]                  // field frequency, 1x1 real
flq.Order=[0]                      // order of photon process, 1x1, int
flq.Amplitude=[1,1,1]              // AC amplitude, 1x1 / 1x2 / 1x3, real
flq.Phase=[0,0,0]                  // AC phase 1x1 / 1x2 / 1x3, real

============= PiLib Variable =============
flq.state_info, @full, [state_label, order, site, identifier, l, SubOrb]
ORDER=    0, SIZE=[   16,    6], TYPE=INTEGER

       1       2       3       4       5       6

       1       0       1       1       2       9
       2       0       1       1       2      10
       3       0       1       1       2      11
       4       0       1       1       2      12
       5       0       1       1       2      13
       6       0       1       1       2      14
       7       0       1       1       2      15
       8       0       1       1       2      16
       9       0       1       1       2      17
      10       0       1       1       2      18
      11       0       2       2       1       3
      12       0       2       2       1       4
      13       0       2       2       1       5
      14       0       2       2       1       6
      15       0       2       2       1       7
      16       0       2       2       1       8

============= PiLib Variable =============
flq.H_onsite(1), @t-sp, Floquet H_onsite of order 0
ORDER=    1, SIZE=[   17,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       1   -0.798910  0.000000
       2       2   -0.798680  0.000000
       3       3   -0.798650  0.000000
       4       4    0.346530  0.000000
       5       5    0.353720  0.000000
       6       6   -0.798910  0.000000
       7       7   -0.798680  0.000000
       8       8   -0.798650  0.000000
       9       9    0.346530  0.000000
      10      10    0.353720  0.000000
      11      11   -0.300000  0.000000
      12      12   -0.300000  0.000000
      13      13   -0.300000  0.000000
      14      14   -0.300000  0.000000
      15      15   -0.300000  0.000000
      16      16   -0.300000  0.000000

============= PiLib Variable =============
flq.hop_size, @full, size of flq.hop_mat [order+1,sublatt,hop_mat_size]
ORDER=    0, SIZE=[    2,    5], TYPE=INTEGER

       1       2       3       4       5

       1       1      16      16       6
       1       2      16      16       6

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,1), @a-sp, Floquet hop_mat(1)(:,:,1) of order 0
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3      13    0.007181  0.000000
       4      11    0.012439  0.000000
       5      11   -0.007181  0.000000
       8      16    0.007181  0.000000
       9      14    0.012439  0.000000
      10      14   -0.007181  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,2), @a-sp, Floquet hop_mat(1)(:,:,2) of order 0
ORDER=    1, SIZE=[    9,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1      11    0.007181  0.000000
       2      13    0.007181  0.000000
       4      12   -0.012439  0.000000
       5      12   -0.007181  0.000000
       6      14    0.007181  0.000000
       7      16    0.007181  0.000000
       9      15   -0.012439  0.000000
      10      15   -0.007181  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,3), @a-sp, Floquet hop_mat(1)(:,:,3) of order 0
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       2      12   -0.007181  0.000000
       3      11   -0.007181  0.000000
       5      13   -0.014363  0.000000
       7      15   -0.007181  0.000000
       8      14   -0.007181  0.000000
      10      16   -0.014363  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,4), @a-sp, Floquet hop_mat(1)(:,:,4) of order 0
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       2      12    0.007181  0.000000
       3      11    0.007181  0.000000
       5      13    0.014363  0.000000
       7      15    0.007181  0.000000
       8      14    0.007181  0.000000
      10      16    0.014363  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,5), @a-sp, Floquet hop_mat(1)(:,:,5) of order 0
ORDER=    1, SIZE=[    9,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1      11   -0.007181  0.000000
       2      13   -0.007181  0.000000
       4      12    0.012439  0.000000
       5      12    0.007181  0.000000
       6      14   -0.007181  0.000000
       7      16   -0.007181  0.000000
       9      15    0.012439  0.000000
      10      15    0.007181  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,6), @a-sp, Floquet hop_mat(1)(:,:,6) of order 0
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3      13   -0.007181  0.000000
       4      11   -0.012439  0.000000
       5      11    0.007181  0.000000
       8      16   -0.007181  0.000000
       9      14   -0.012439  0.000000
      10      14    0.007181  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,1), @a-sp, Floquet hop_mat(2)(:,:,1) of order 0
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       4   -0.012439  0.000000
      11       5    0.007181  0.000000
      13       3   -0.007181  0.000000
      14       9   -0.012439  0.000000
      14      10    0.007181  0.000000
      16       8   -0.007181  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,2), @a-sp, Floquet hop_mat(2)(:,:,2) of order 0
ORDER=    1, SIZE=[    9,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       1   -0.007181  0.000000
      12       4    0.012439  0.000000
      12       5    0.007181  0.000000
      13       2   -0.007181  0.000000
      14       6   -0.007181  0.000000
      15       9    0.012439  0.000000
      15      10    0.007181  0.000000
      16       7   -0.007181  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,3), @a-sp, Floquet hop_mat(2)(:,:,3) of order 0
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       3    0.007181  0.000000
      12       2    0.007181  0.000000
      13       5    0.014363  0.000000
      14       8    0.007181  0.000000
      15       7    0.007181  0.000000
      16      10    0.014363  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,4), @a-sp, Floquet hop_mat(2)(:,:,4) of order 0
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       3   -0.007181  0.000000
      12       2   -0.007181  0.000000
      13       5   -0.014363  0.000000
      14       8   -0.007181  0.000000
      15       7   -0.007181  0.000000
      16      10   -0.014363  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,5), @a-sp, Floquet hop_mat(2)(:,:,5) of order 0
ORDER=    1, SIZE=[    9,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       1    0.007181  0.000000
      12       4   -0.012439  0.000000
      12       5   -0.007181  0.000000
      13       2    0.007181  0.000000
      14       6    0.007181  0.000000
      15       9   -0.012439  0.000000
      15      10   -0.007181  0.000000
      16       7    0.007181  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,6), @a-sp, Floquet hop_mat(2)(:,:,6) of order 0
ORDER=    1, SIZE=[    7,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       4    0.012439  0.000000
      11       5   -0.007181  0.000000
      13       3    0.007181  0.000000
      14       9    0.012439  0.000000
      14      10   -0.007181  0.000000
      16       8    0.007181  0.000000
