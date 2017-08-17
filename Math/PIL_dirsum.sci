// **** Purpose ****
// calculate the direct sum of two matrix
// **** Variables ****
// A: nxn, real or complex
// <= matrix A
// B: nxn, real or complex
// <= matrix B
// C: nxn, real or complex
// => direct sum of A and B
// **** Version ****
// 05/01/2014
// **** Comment ****

function [C]=PIL_dirsum(A,B)
    A_c=length(A(1,:));
    A_r=length(A(:,1));
    B_c=length(B(1,:));
    B_r=length(B(:,1));
    C=zeros(A_r+B_r,A_c+B_c);
    C(1:A_r,1:A_c)=A;
    C(A_r+1:A_r+B_r,A_c+1:A_c+B_c)=B;
endfunction
