// **** Purpose ****
// try to block diagonal a input matrix
// **** Variables ****
// [M_in]: n x n, complex
// <= the input matrix
// [M_group]: n x 2, integer
// => output [group label, state label]
// [M_out]: nxn, complex
// => output matrix
// **** Version ****
// 05/01/2014
// **** Comment ****

function [M_group,M_out]=PIL_blkdiag(M_in)
    M_sp=PIL_sparse(M_in,'full','all');
    // find group & relabel
    M_group=[1:M_sp(1,1);1:M_sp(1,1)]';
    for n=2:length(M_sp(:,1))
        M_group(find(M_group(:,1)==M_sp(n,1)),1)=M_sp(n,2);
    end
    M_group=gsort(M_group,'lr','i');
    
    // generate unitary transformation
    U=zeros(M_sp(1),M_sp(2));
    for n=1:M_sp(1)
        U(M_group(n,2),n)=1;
    end
    M_out=U'*M_in*U
endfunction

// example:
// M_in=[1 0 1 0 1; 0 1 0 1 0; 1 0 1 0 1; 0 1 0 1 0; 1 0 1 0 1];
// M_out  =
// 
//    1.    1.    1.    0.    0.  
//    1.    1.    1.    0.    0.  
//    1.    1.    1.    0.    0.  
//    0.    0.    0.    1.    1.  
//    0.    0.    0.    1.    1.  
