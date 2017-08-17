// **** Purpose ****
// calculate the commutation relation 
// **** Variables ****
// A: nxn, real or complex
// <= matrix A
// B: nxn, real or complex
// <= matrix B
// cal_type: 1x1, char, default='c'
// <= 'a' (anticommuntation), 'c'(communtation)
// C: nxn, real or complex
// => communtation matrix
// **** Version ****
// 05/01/2014
// **** Comment ****

function [C]=PIL_commutator(A,B,cal_type)
    [lhs,rhs]=argn();
    if rhs==2 then
        cal_type='c'
    end
    select cal_type
    case 'a'
        C=A*B+B*A;
    case 'c'
        C=A*B-B*A;
    end
endfunction
