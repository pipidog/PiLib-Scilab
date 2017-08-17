// **** Purpose ****
// Recover matrix from coeffiients of a tensor operator basis 
// **** Variables ****
// T_coff: n^2 x 1, real or complex 
// <= the expansion coefficients of a tensor opeator basis 
// basis_op: n x n x n^2, real or complex 
// <= the tensor opeatro basis, a three indices variable 
// M_out: nxn, real or complex
// => the recombined matrix
// **** Version ****
// 05/01/2014
// **** Comment ****
// the function name "TO=tensor operator"

function [M_out]=PIL_TO_recomb(T_coff,basis_op)
    J=(sqrt(length(T_coff))-1)/2;
    [lhs,rhs]=argn()
    
    if rhs==1
        basis_op='s';
    end
    
    if type(basis_op)==10
        T=mylib_tensor_op(J,basis_op);
    else
        T=basis_op;
    end
    T_size=size(T);
    M_out=zeros(T_size(1),T_size(2));
    for n=1:T_size(3)
        M_out=M_out+T_coff(n)*T(:,:,n);
    end
endfunction
