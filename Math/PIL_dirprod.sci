// **** Purpose ****
// calculate the direct product
// **** Variables ****
// A: nxn, real or complex
// <= matrix A
// B: nxn, real or complex
// <= matrix B
// C: nxn, real or complex
// => A direct prod of B
// **** Version ****
// 05/01/2014
// **** Comment ****

function [C]=PIL_dirprod(A,B)
 A_dim=size(A);
 B_dim=size(B);
 C=zeros(A_dim(1)*B_dim(1),A_dim(2)*B_dim(2));
 for n=1:A_dim(1)
     for m=1:A_dim(2)
         C((n-1)*B_dim(1)+1:n*B_dim(1),(m-1)*B_dim(2)+1:m*B_dim(2))=A(n,m)*B;
     end
 end
endfunction
