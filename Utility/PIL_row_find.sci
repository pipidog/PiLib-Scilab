// **** Purpose ****
// This function finds whether there is a row of matrix A matches B
// **** Variables ****
// [A]: complex, NxM 
// <= full NxM numerical matrix 1  
// [B]: complex, 1xM 
// <= single row 1xM numerical matrix 2
// [find_index]: index of the row in A that matches B
// **** Version ****
// 11/11/2014
function find_index=PIL_row_find(A,B,criterion)
    A_size=size(A); B_size=size(B)
    if (A_size(2)~=B_size(2)) | (B_size(1)~=1) then
        disp('Error: PIL_row_find, A and B dimension inconsistent!');
        abort
    end
    [lhs,rhs]=argn()
    if rhs==2 then
        criterion=1e-5;
    end
    C=clean(abs(A-repmat(B,length(A(:,1)),1)),criterion);
    D=zeros(A_size(1),1);
    for n=1:A_size(2)
        D=D+C(:,n);
    end
    D=clean(D);
    find_index=find(D<=1e-7);
endfunction
