// **** Purpose ****
// calculate the cross product of two vectors
// **** Variables ****
// V1: nx1, real
// <= vector 1
// V2: nx1, real
// <= vector 2
// V3: nx1, real
// => cross product of V1xV2
// **** Version ****
// 05/01/2014
// **** Comment ****

function V3=PIL_crossprod(V1, V2)
 V3 = [V1(2) * V2(3) - V1(3) * V2(2) ; 
        V1(3) * V2(1) - V1(1) * V2(3)
        V1(1) * V2(2) - V1(2) * V2(1)];
endfunction
