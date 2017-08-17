flq.Frequency=[10]                  // field frequency, 1x1 real
flq.Order=[1]                      // order of photon process, 1x1, int
flq.Amplitude=[1,0]                // AC amplitude, 1x1 / 1x2 / 1x3, real
flq.Phase=[%pi/2,0]                // AC phase 1x1 / 1x2 / 1x3, real

============= PiLib Variable =============
flq.state_info, @full, [state_label, order, site, identifier, l, SubOrb]
ORDER=    0, SIZE=[   48,    6], TYPE=INTEGER

       1       2       3       4       5       6

       1      -1       1       1       0       1
       2      -1       1       1       0       2
       3      -1       2       2       0       1
       4      -1       2       2       0       2
       5      -1       3       3       0       1
       6      -1       3       3       0       2
       7      -1       4       4       0       1
       8      -1       4       4       0       2
       9      -1       5       5       0       1
      10      -1       5       5       0       2
      11      -1       6       6       0       1
      12      -1       6       6       0       2
      13      -1       7       7       0       1
      14      -1       7       7       0       2
      15      -1       8       8       0       1
      16      -1       8       8       0       2
      17       0       1       1       0       1
      18       0       1       1       0       2
      19       0       2       2       0       1
      20       0       2       2       0       2
      21       0       3       3       0       1
      22       0       3       3       0       2
      23       0       4       4       0       1
      24       0       4       4       0       2
      25       0       5       5       0       1
      26       0       5       5       0       2
      27       0       6       6       0       1
      28       0       6       6       0       2
      29       0       7       7       0       1
      30       0       7       7       0       2
      31       0       8       8       0       1
      32       0       8       8       0       2
      33       1       1       1       0       1
      34       1       1       1       0       2
      35       1       2       2       0       1
      36       1       2       2       0       2
      37       1       3       3       0       1
      38       1       3       3       0       2
      39       1       4       4       0       1
      40       1       4       4       0       2
      41       1       5       5       0       1
      42       1       5       5       0       2
      43       1       6       6       0       1
      44       1       6       6       0       2
      45       1       7       7       0       1
      46       1       7       7       0       2
      47       1       8       8       0       1
      48       1       8       8       0       2

============= PiLib Variable =============
flq.H_onsite(1), @t-sp, Floquet H_onsite of order 0
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
flq.H_onsite(2), @t-sp, Floquet H_onsite of order 1
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
flq.hop_size, @full, size of flq.hop_mat [order+1,sublatt,hop_mat_size]
ORDER=    0, SIZE=[   16,    5], TYPE=INTEGER

       1       2       3       4       5

       1       1      16      16       2
       1       2      16      16       3
       1       3      16      16       3
       1       4      16      16       3
       1       5      16      16       3
       1       6      16      16       3
       1       7      16      16       3
       1       8      16      16       2
       2       1      16      16       2
       2       2      16      16       3
       2       3      16      16       3
       2       4      16      16       3
       2       5      16      16       3
       2       6      16      16       3
       2       7      16      16       3
       2       8      16      16       2

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,1), @a-sp, Floquet hop_mat(1)(:,:,1) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       3    0.082111  0.000000
       2       4    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,2), @a-sp, Floquet hop_mat(1)(:,:,2) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       3    0.082111  0.000000
       2       4    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,1), @a-sp, Floquet hop_mat(2)(:,:,1) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3       1    0.082111  0.000000
       4       2    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,2), @a-sp, Floquet hop_mat(2)(:,:,2) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3       1    0.082111  0.000000
       4       2    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,3), @a-sp, Floquet hop_mat(2)(:,:,3) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3       5    0.100000  0.000000
       4       6    0.100000  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(3)(:,:,1), @a-sp, Floquet hop_mat(3)(:,:,1) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       5       3    0.100000  0.000000
       6       4    0.100000  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(3)(:,:,2), @a-sp, Floquet hop_mat(3)(:,:,2) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       5       7    0.082111  0.000000
       6       8    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(3)(:,:,3), @a-sp, Floquet hop_mat(3)(:,:,3) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       5       7    0.082111  0.000000
       6       8    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(4)(:,:,1), @a-sp, Floquet hop_mat(4)(:,:,1) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       7       5    0.082111  0.000000
       8       6    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(4)(:,:,2), @a-sp, Floquet hop_mat(4)(:,:,2) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       7       5    0.082111  0.000000
       8       6    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(4)(:,:,3), @a-sp, Floquet hop_mat(4)(:,:,3) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       7       9    0.100000  0.000000
       8      10    0.100000  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(5)(:,:,1), @a-sp, Floquet hop_mat(5)(:,:,1) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       9       7    0.100000  0.000000
      10       8    0.100000  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(5)(:,:,2), @a-sp, Floquet hop_mat(5)(:,:,2) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       9      11    0.082111  0.000000
      10      12    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(5)(:,:,3), @a-sp, Floquet hop_mat(5)(:,:,3) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       9      11    0.082111  0.000000
      10      12    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(6)(:,:,1), @a-sp, Floquet hop_mat(6)(:,:,1) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       9    0.082111  0.000000
      12      10    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(6)(:,:,2), @a-sp, Floquet hop_mat(6)(:,:,2) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       9    0.082111  0.000000
      12      10    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(6)(:,:,3), @a-sp, Floquet hop_mat(6)(:,:,3) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11      13    0.100000  0.000000
      12      14    0.100000  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(7)(:,:,1), @a-sp, Floquet hop_mat(7)(:,:,1) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      13      11    0.100000  0.000000
      14      12    0.100000  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(7)(:,:,2), @a-sp, Floquet hop_mat(7)(:,:,2) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      13      15    0.082111  0.000000
      14      16    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(7)(:,:,3), @a-sp, Floquet hop_mat(7)(:,:,3) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      13      15    0.082111  0.000000
      14      16    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(8)(:,:,1), @a-sp, Floquet hop_mat(8)(:,:,1) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      15      13    0.082111  0.000000
      16      14    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(8)(:,:,2), @a-sp, Floquet hop_mat(8)(:,:,2) of order 0
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      15      13    0.082111  0.000000
      16      14    0.082111  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,1), @a-sp, Floquet hop_mat(1)(:,:,1) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       3   -0.000000 -0.039367
       2       4   -0.000000 -0.039367

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,2), @a-sp, Floquet hop_mat(1)(:,:,2) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       1       3    0.000000  0.039367
       2       4    0.000000  0.039367

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,1), @a-sp, Floquet hop_mat(2)(:,:,1) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3       1   -0.000000 -0.039367
       4       2   -0.000000 -0.039367

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,2), @a-sp, Floquet hop_mat(2)(:,:,2) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       3       1    0.000000  0.039367
       4       2    0.000000  0.039367

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,3), @a-sp, Floquet hop_mat(2)(:,:,3) of order 1
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(3)(:,:,1), @a-sp, Floquet hop_mat(3)(:,:,1) of order 1
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(3)(:,:,2), @a-sp, Floquet hop_mat(3)(:,:,2) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       5       7   -0.000000 -0.039367
       6       8   -0.000000 -0.039367

============= PiLib Variable =============
flq.hop_mat(2)(3)(:,:,3), @a-sp, Floquet hop_mat(3)(:,:,3) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       5       7    0.000000  0.039367
       6       8    0.000000  0.039367

============= PiLib Variable =============
flq.hop_mat(2)(4)(:,:,1), @a-sp, Floquet hop_mat(4)(:,:,1) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       7       5   -0.000000 -0.039367
       8       6   -0.000000 -0.039367

============= PiLib Variable =============
flq.hop_mat(2)(4)(:,:,2), @a-sp, Floquet hop_mat(4)(:,:,2) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       7       5    0.000000  0.039367
       8       6    0.000000  0.039367

============= PiLib Variable =============
flq.hop_mat(2)(4)(:,:,3), @a-sp, Floquet hop_mat(4)(:,:,3) of order 1
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(5)(:,:,1), @a-sp, Floquet hop_mat(5)(:,:,1) of order 1
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(5)(:,:,2), @a-sp, Floquet hop_mat(5)(:,:,2) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       9      11   -0.000000 -0.039367
      10      12   -0.000000 -0.039367

============= PiLib Variable =============
flq.hop_mat(2)(5)(:,:,3), @a-sp, Floquet hop_mat(5)(:,:,3) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
       9      11    0.000000  0.039367
      10      12    0.000000  0.039367

============= PiLib Variable =============
flq.hop_mat(2)(6)(:,:,1), @a-sp, Floquet hop_mat(6)(:,:,1) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       9   -0.000000 -0.039367
      12      10   -0.000000 -0.039367

============= PiLib Variable =============
flq.hop_mat(2)(6)(:,:,2), @a-sp, Floquet hop_mat(6)(:,:,2) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      11       9    0.000000  0.039367
      12      10    0.000000  0.039367

============= PiLib Variable =============
flq.hop_mat(2)(6)(:,:,3), @a-sp, Floquet hop_mat(6)(:,:,3) of order 1
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(7)(:,:,1), @a-sp, Floquet hop_mat(7)(:,:,1) of order 1
ORDER=    1, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(7)(:,:,2), @a-sp, Floquet hop_mat(7)(:,:,2) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      13      15   -0.000000 -0.039367
      14      16   -0.000000 -0.039367

============= PiLib Variable =============
flq.hop_mat(2)(7)(:,:,3), @a-sp, Floquet hop_mat(7)(:,:,3) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      13      15    0.000000  0.039367
      14      16    0.000000  0.039367

============= PiLib Variable =============
flq.hop_mat(2)(8)(:,:,1), @a-sp, Floquet hop_mat(8)(:,:,1) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      15      13   -0.000000 -0.039367
      16      14   -0.000000 -0.039367

============= PiLib Variable =============
flq.hop_mat(2)(8)(:,:,2), @a-sp, Floquet hop_mat(8)(:,:,2) of order 1
ORDER=    1, SIZE=[    3,    3], TYPE=SPARSE

       1       2                     3

      16      16    0.000000  0.000000
      15      13    0.000000  0.039367
      16      14    0.000000  0.039367
