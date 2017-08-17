// **** Purpose ****
// This function split a word into separate strings by identifying '#' 
// **** Variables ****
// [A]: string
// <= a line of strings
// [B]: string
// => separated strings
// **** Version ****
// 06/07/2014 1st version
// **** Comment ****
// 1. This is a PiLib IO function. In PiLib, any string ending with '#' is
//    idenfitied as a string element in a string matrix. 
// 2. The leading and trailing blanks (and tabs) of strings will be automatically
//    trimmed.

function B=PIL_str_split(A)
    A_len=length(A);
    A_char=strsplit(A);
    A_pound=find(A_char=='#');
    // generate string range
    A_range=zeros(length(A_pound),2);
    A_range(1,1)=1
    A_range($,2)=A_pound($)-1;
    
    for n=1:length(A_pound)-1
        A_range(n,2)=A_pound(n)-1;
        A_range(n+1,1)=A_pound(n)+1
    end
    B=emptystr(1,length(A_pound));
    for n=1:length(A_pound)
        for m=A_range(n,1):A_range(n,2)
            B(n)=B(n)+A_char(m);
        end
        B(n)=stripblanks(B(n));
    end
endfunction
