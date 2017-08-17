// **** Purpose ****
// calculate the Kroneker delta 
// **** Variables ****
// n: 1x1, integer
// <= value n
// n: 1x1, integer
// <= value n
// p: 1x1, integer, 0 or 1
// => Kronker delta of a,b
// **** Version ****
// 05/01/2014
// **** Comment ****

function p=PIL_kdelta(n,m)
    if n==m then
        p=1;
    else
        p=0;
    end
endfunction
