// **** Purpose ****
// expand a vector using another basis set
// **** Variables ****
// V: nx1, real or complex
// <= vector to expan
// M: nxn, real or complex
// <= the basis set matrix
// Expan_coef: nx1, integer
// => the expansion coefficients
// **** Version ****
// 05/01/2014
// **** Comment ****
// 1. V can be row or column vectros. However, M must be defined by 
//    column vectros!
// 2. M basis set can be non-orthnormal

function Expan_coef=PIL_linexpan(V,M)
    V_size=size(V);
    if V_size(2)~=1 then
        V=V';
    end
    Expan_coef=inv(M)*V;
endfunction
