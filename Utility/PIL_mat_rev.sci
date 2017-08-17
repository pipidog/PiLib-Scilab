// **** Purpose ****
// reverse the order of states that expand a square matrix 
// **** Variables ****
// [M]: n x n, complex
// <= the input matrix
// [M_rev]: n x x, complex
// => the reversed matrix
// **** Version ****
// 05/01/2014
// **** Comment ****

function M_rev=PIL_mat_rev(M)
    M_rev=flipdim(flipdim(M,1),2);
endfunction
