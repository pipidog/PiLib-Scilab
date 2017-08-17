// **** Purpose ****
// calculate the trace of a matrix
// **** Variables ****
// A: nxn, complex or real
// <= input matrix
// output: 1x1, complex or real
// => trace of A 
// **** Version ****
// 05/01/2014
// **** Comment ****

function output=PIL_trace(A)
    output=sum(diag(A));
endfunction
