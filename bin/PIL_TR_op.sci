// **** Purpose ****
// This function perform the Time-Reveresal operation of a given state
// in PiLab format. Input can be normal states or floquet states. 
// **** Variables ****
// [A]: Nx1 / NxM, 
// <= eigenvector or eigenmatrix
// => TR operated eigenvector or eigenmatrix
// [TR_pair]: variables generated by PIL_TR_pair
// **** Version ****
// 04/30/2015 first built
// **** Comment ****
// 1. This function perform TR operation on flq or hop states. 
// The operation requires hop.Basis='c' and the flq.Phase=[0,0]
// You must check it before using this code ! 
// 2. The variable TR_pair tells this code the does the state paired with
// their time-reversal partner. It should be obtaine from PIL_TR_pair  
function [A]=PIL_TR_op(A,TR_pair)
    // check TR_pair
    if length(TR_pair(:,1))~=length(A(:,1))/2 then
        disp('Error: PIL_TR_op, size of TR_pair is inconsistent with A');
        abort;
    end
    // check input object A
    A_size=size(A);
    if A_size(1)==1
        disp('Warning: PIL_TR_op, A is converted to coulmn vector!');
        A=A.';
    end
    // TR operation
    for n=1:length(TR_pair(:,1))
        A([TR_pair(n,:)],:)=conj(A(flipdim([TR_pair(n,:)],2),:));
        A(TR_pair(n,2),:)=-A(TR_pair(n,2),:);
    end
endfunction


