// **** Purpose ****
// Convert the input sublattice from reduced coordinates to cartesian
// coordinate and vice versa. 
// **** Variables ****
// [primitive]: 3x3, real
// <= primitive row vectors 
// [sublat_in]: nx3, real
// <= input coordinates of sublatice
// [inp_type]: 1x1, string, 'red' / 'cart' / 'cart-shift'
// <= specify type of input sublat, reduced or cartesian coordinate
//    'cart-shift' will force the output reduced coordinate in the
//    center [0,0,0] unitcell. 
// [sublat_out]: nx3, real
// => output coordinates of sublatice
// **** Version ****
// 01/07/2016
// **** Comment ****
// reduced coordniate: c1*a1+c2*a2+c3*a3 =>[c1,c2,c3]
// cartesian coordinate: [x,y,z]

function sublat_out=PIL_red_cart_conv(primitive, sublat_in, inp_type)
    tot_at=length(sublat_in(:,1));
    select inp_type
    case 'red'
        sublat_out=sublat_in*primitive
    case 'cart'
        for n=1:tot_at
            sublat_out(n,:)=(PIL_linexpan(sublat_in(n,:),primitive'))';
        end
    case 'cart-shift'
        for n=1:tot_at
            sublat_out(n,:)=(PIL_linexpan(sublat_in(n,:),primitive'))';
        end
        sublat_out=sublat_out-fix(sublat_out);
    end
endfunction
