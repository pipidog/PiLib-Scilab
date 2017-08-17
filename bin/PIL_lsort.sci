// **** Purpose ****
// lexicographic row/column sort of a input matrix based on the 
// assigned row/column sequence.
// **** Variables ****
// [A_in]: complex, NxM 
// <= numerical matrix 1   
// [r_c_type]: 1x1, string, 'r' / 'c'
// <= 'r': row, 'c': column
// [data_seq]: 1xN / 1xM, int
// <= the sequence of the sort 
// [direction]: 1x1, string: 'i' / 'd'
// <= 'i' -> increasing, 'd' -> decending
// [A_out]: complex, NxM
// => sorted Matrix
// **** Version ****
// 5/19/2015: first built
// **** Comment ****
function A_out=PIL_lsort(A_in,r_c_type,data_seq,direction)
    [lhs,rhs]=argn()
    if rhs==3
        direction='i'
    end
    
    data_seq=cat(1,data_seq,[1:length(data_seq)]);
    select r_c_type 
    case 'c'
        if length(data_seq(1,:))~=length(A_in(1,:))
            disp('Error: PIL_lsort, size of data_seq ~= # of columns in A_in')
            abort;
        end  
        A_out=A_in(:,[data_seq(1,:)]);
        A_out=gsort(A_out,'lr',direction);
        data_seq=gsort(data_seq,'lc','i');
        A_out=A_out(:,[data_seq(2,:)]);
    case 'r'
        if length(data_seq(1,:))~=length(A_in(:,1))
            disp('Error: PIL_lsort, size of data_seq ~= # of rows in A_in')
            abort;
        end
        A_out=A_in([data_seq(1,:)],:);
        A_out=gsort(A_out,'lc',direction);
        data_seq=gsort(data_seq,'lc','i');
        A_out=A_out([data_seq(2,:)],:);
    end
endfunction
