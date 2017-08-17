// **** Purpose ****
// calculate the 3j symbol
// **** Variables ****
// [j1]: 1x1, iteger or half-integer
// <= j1 parameter
// [m1]: 1x1, iteger or half-integer
// [val]: 1x1, real or complex
// => value of 3j symbol 
// => cobinatorial
// **** Version ****
// 05/01/2014
// **** Comment ****
// This program provides you the values of 3-j symbol
// see "Wigner 3-j symbols" (wiki)

function val=PIL_tri_j_sym(j1,j2,j3,m1,m2,m3)
    val=( ((-1)^(j1-j2-m3))/sqrt(2*j3+1) )*PIL_cg_coff(j1,j2,m1,m2,j3,-m3)
endfunction
