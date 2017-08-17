// **** Purpose ****
// calculates the voulmn of a unit cell for a given primitive cell
// **** Variables ****
// [primitive_cell]: 3x3, real
// <= primitice vectros in row
// [V]: 1x1, real
// => the volumn
// **** Version ****
// 05/01/2014
// **** Comment ****
function V=PIL_V_primitive(primitive_cell)
    a_r=primitive_cell;
    V=abs(a_r(1,:)*PIL_crossprod(a_r(2,:)',a_r(3,:)'));
endfunction
