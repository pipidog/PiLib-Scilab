// **** Purpose ****
// This code calculates the Gaunt cofficients
// **** Variables ****
// [j1]: 1x1, integer or half-integer
// <= j1 parameter
// [m1]: 1x1, integer or half-integer
// <= m1 parameter
// [G_coff]: 1x1, real or complex
// => Gaunt co
// **** Version ****
// Mar/13/2014 First Built
// **** Comment ****
// Ref: http://theoretical-physics.net/dev/src/math/spherical-harmonics.html
function G_coff=PIL_Gaunt_coff(j1,j2,j3,m1,m2,m3)
    G_coff=(-1)^(-m1)*sqrt((2*j1+1)*(2*j2+1)*(2*j3+1)/(4*%pi))...
    *PIL_tri_j_sym(j1,j2,j3,0,0,0)*PIL_tri_j_sym(j1,j2,j3,-m1,m2,m3);
endfunction















