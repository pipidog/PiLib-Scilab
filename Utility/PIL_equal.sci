// **** Purpose ****
// This function identifies whether two input numerical objects equal
// **** Variables ****
// [A]: complex, NxM 
// <= numerical matrix 1  
// [B]: complex, NxM 
// <= numerical matrix 2
// [BL_val]: boolean value of their comparsion
// **** Version ****
// 11/11/2014
function BL_val=PIL_equal(A,B,criterion)
    [lhs,rhs]=argn()
    if rhs==2 then
        criterion=1e-6;
    end
    if sum(abs(size(A)-size(B)))~=0 then
        disp('Error: PIL_equal, input objects inconsistent!');
        abort
    end
    BL_val=max(abs(A-B)) < criterion;
endfunction
