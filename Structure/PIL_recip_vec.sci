// **** Purpose ****
// This program gives you the primitive vectors of reciprocal lattice.
// **** Variables ****
// [primitive_cell]: 3x3, real
// <= primitive cell in row vectors
// [b]: 3x3, real
// => reciprocal lattice in row vectros
// **** Version ****
// 05/01/2014 first version
// 06/03/2014 fully rewrite to accept all dimension
// **** Comment ****

function b=PIL_recip_vec(primitive_cell)
    select length(primitive_cell(:,1))
    case 3
        b1=2*%pi*inv(primitive_cell)*[1;0;0];
        b2=2*%pi*inv(primitive_cell)*[0;1;0];
        b3=2*%pi*inv(primitive_cell)*[0;0;1];
        b=[b1';b2';b3'];
    case 2
        b1=2*%pi*inv(primitive_cell)*[1;0];
        b2=2*%pi*inv(primitive_cell)*[0;1];
        b=[b1';b2'];
    case 1
        b1=2*%pi*inv(primitive_cell)*[1];
        b=[b1'];
    end
endfunction

