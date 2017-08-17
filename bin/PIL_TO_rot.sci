// **** Purpose ****
// rotate an input tensor operator to different quantization axis 
// **** Variables ****
// T_in: n x n x n^2, real or complex
// <= input tensor operator basis
// Euler_angles: 1x3, real
// <= Euler angle (alpha, beta, gamma ; see Tzeng V-II p353 )
// T_out: n x n x n^2, real or complex
// => output tensor opeartor basis
// **** Version ****
// 05/01/2014
// **** Comment ****
// 1. This function rotates all tensor operators with an specific Euler angles
// 2. Notice: all your tensor operators should be presented in |J,MJ> basis.
// 3. T-operators whose quantization axis are defined in r, then this function
//    rotate its quantization axis to r'=D*r.
// 4. T_in is the whole basis

function [T_out]=PIL_TO_rot(T_in,Euler_angles)
    J=((length(T_in(:,1,1)))-1)/2;
    D=PIL_D_mat(J,Euler_angles);
    T_size=size(T_in);
    for n=1:T_size(3)
        T_out(:,:,n)=D*T_in(:,:,n)*D'
    end
endfunction
