// **** Purpose ****
// generates the rotation matrix of a classical system. The rotaton is 
// defined by rotate along an axis with angle phi rather than Euler angles 
// **** Variables ****
// [u]: 1x3, real
// <= the axis vector of the rotation. only direction matters. the code 
//    will automatically normalize its length in calculations.
// [phi]: 1x1, real
// <= the rotation angles between +Pi ~ -%Pi
// [M]: 3x3, real
// => The corresponding rotation matrix
// **** Version ****
// 01/21/2014
// **** Comment ****
// A few useful applications of this approach:
// 1). For two unit vectors A, B, if we want to generate a rotation matrix
//     that rotates A to B, i.e. MA=B. then u and phi are just:
//     u=(A+B)/2, phi=%pi
//     (Note: A and B must be unit vectors!!)
// 2). If one wants to specify this rotation matrix in term of Euler
//     angles, pass the result to PIL_Euler_finder. 

function M=PIL_rot_axis(u,phi)
    u=u/norm(u);
    M=[cos(phi)+u(1)^2*(1-cos(phi)),  u(1)*u(2)*(1-cos(phi))-u(3)*sin(phi),  u(1)*u(3)*(1-cos(phi))+u(2)*sin(phi);
    u(2)*u(1)*(1-cos(phi))+u(3)*sin(phi), cos(phi)+u(2)^2*(1-cos(phi)),  u(2)*u(3)*(1-cos(phi))-u(1)*sin(phi);
    u(3)*u(1)*(1-cos(phi))-u(2)*sin(phi),  u(3)*u(2)*(1-cos(phi))+u(1)*sin(phi),  cos(phi)+u(3)^2*(1-cos(phi))]
endfunction
