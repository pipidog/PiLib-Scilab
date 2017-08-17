flq.Frequency=[10]                  // field frequency, 1x1 real
flq.Order=[1]                      // order of photon process, 1x1, int
flq.Amplitude=[1,1]                // AC amplitude, 1x1 / 1x2 / 1x3, real
flq.Phase=[0,%pi/2]                    // AC phase 1x1 / 1x2 / 1x3, real

============= PiLib Variable =============
flq.state_info, @full, [state_label, order, site, identifier, l, SubOrb]
ORDER=    0, SIZE=[    6,    6], TYPE=INTEGER

       1       2       3       4       5       6

       1      -1       1       1       1       5
       2      -1       2       2       1       5
       3       0       1       1       1       5
       4       0       2       2       1       5
       5       1       1       1       1       5
       6       1       2       2       1       5

============= PiLib Variable =============
flq.H_onsite(1), @t-sp, Floquet H_onsite of order 0
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000

============= PiLib Variable =============
flq.H_onsite(2), @t-sp, Floquet H_onsite of order 1
ORDER=    0, SIZE=[    1,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000

============= PiLib Variable =============
flq.hop_size, @full, size of flq.hop_mat [order+1,sublatt,hop_mat_size]
ORDER=    0, SIZE=[    4,    5], TYPE=INTEGER

       1       2       3       4       5

       1       1       2       2       9
       1       2       2       2       9
       2       1       2       2       9
       2       2       2       2       9

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,1), @a-sp, Floquet hop_mat(1)(:,:,1) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       2    0.765198  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,2), @a-sp, Floquet hop_mat(1)(:,:,2) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       2    0.765198  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,3), @a-sp, Floquet hop_mat(1)(:,:,3) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       2    0.765198  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,4), @a-sp, Floquet hop_mat(1)(:,:,4) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000  0.037944

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,5), @a-sp, Floquet hop_mat(1)(:,:,5) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000 -0.037944

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,6), @a-sp, Floquet hop_mat(1)(:,:,6) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000 -0.037944

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,7), @a-sp, Floquet hop_mat(1)(:,:,7) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000  0.037944

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,8), @a-sp, Floquet hop_mat(1)(:,:,8) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000  0.037944

============= PiLib Variable =============
flq.hop_mat(1)(1)(:,:,9), @a-sp, Floquet hop_mat(1)(:,:,9) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.000000 -0.037944

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,1), @a-sp, Floquet hop_mat(2)(:,:,1) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       1    0.765198  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,2), @a-sp, Floquet hop_mat(2)(:,:,2) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       1    0.765198  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,3), @a-sp, Floquet hop_mat(2)(:,:,3) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       1    0.765198  0.000000

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,4), @a-sp, Floquet hop_mat(2)(:,:,4) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000 -0.037944

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,5), @a-sp, Floquet hop_mat(2)(:,:,5) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000  0.037944

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,6), @a-sp, Floquet hop_mat(2)(:,:,6) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000  0.037944

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,7), @a-sp, Floquet hop_mat(2)(:,:,7) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000 -0.037944

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,8), @a-sp, Floquet hop_mat(2)(:,:,8) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000 -0.037944

============= PiLib Variable =============
flq.hop_mat(1)(2)(:,:,9), @a-sp, Floquet hop_mat(2)(:,:,9) of order 0
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.000000  0.037944

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,1), @a-sp, Floquet hop_mat(1)(:,:,1) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       2   -0.220025 -0.381095

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,2), @a-sp, Floquet hop_mat(1)(:,:,2) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       2    0.440051  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,3), @a-sp, Floquet hop_mat(1)(:,:,3) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       2   -0.220025  0.381095

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,4), @a-sp, Floquet hop_mat(1)(:,:,4) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.057942 -0.000000

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,5), @a-sp, Floquet hop_mat(1)(:,:,5) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1   -0.028971 -0.050179

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,6), @a-sp, Floquet hop_mat(1)(:,:,6) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1   -0.028971  0.050179

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,7), @a-sp, Floquet hop_mat(1)(:,:,7) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1   -0.028971  0.050179

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,8), @a-sp, Floquet hop_mat(1)(:,:,8) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1   -0.028971 -0.050179

============= PiLib Variable =============
flq.hop_mat(2)(1)(:,:,9), @a-sp, Floquet hop_mat(1)(:,:,9) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       1       1    0.057942 -0.000000

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,1), @a-sp, Floquet hop_mat(2)(:,:,1) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       1    0.220025 -0.381095

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,2), @a-sp, Floquet hop_mat(2)(:,:,2) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       1   -0.440051  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,3), @a-sp, Floquet hop_mat(2)(:,:,3) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       1    0.220025  0.381095

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,4), @a-sp, Floquet hop_mat(2)(:,:,4) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2   -0.057942  0.000000

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,5), @a-sp, Floquet hop_mat(2)(:,:,5) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.028971  0.050179

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,6), @a-sp, Floquet hop_mat(2)(:,:,6) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.028971 -0.050179

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,7), @a-sp, Floquet hop_mat(2)(:,:,7) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.028971 -0.050179

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,8), @a-sp, Floquet hop_mat(2)(:,:,8) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2    0.028971  0.050179

============= PiLib Variable =============
flq.hop_mat(2)(2)(:,:,9), @a-sp, Floquet hop_mat(2)(:,:,9) of order 1
ORDER=    0, SIZE=[    2,    3], TYPE=SPARSE

       1       2                     3

       2       2    0.000000  0.000000
       2       2   -0.057942  0.000000
