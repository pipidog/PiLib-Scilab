// **** Purpose ****
// Convert the input k-space coordinate from reduced coordinates to 
// cartesian coordinate and vice versa. 
// **** Variables ****
// [prim_vec]: 3x3, real
// <= prim_vec row vectors 
// [kpt_in]: nx3, real
// <= input of the k-point, can be many, in row form
// [inp_type]: 1x1, string, 'red' / 'cart'
// <= specify type of input sublatt, reduced or cartesian coordinate
// [kpt_out]: nx3, real
// => output of the k-point
// **** Version ****
// 01/14/2016
// **** Comment ****
// reduced coordniate: c1*b1+c2*b2+c3*b3 =>[c1,c2,c3]
// cartesian coordinate: [kx,ky,kz]

function kpt_out=PIL_kpt_conv(prim_vec,kpt_in,inp_type)
    tot_k=length(kpt_in(:,1));
    recip_vec=PIL_recip_vec(prim_vec);
    kpt_out=zeros(kpt_in);
    select inp_type
    case 'red'
        for n=1:tot_k // total k-points
            for m=1:length(recip_vec(1,:)) // dimension
                kpt_out(n,:)=kpt_out(n,:)+kpt_in(n,m)*recip_vec(m,:);
            end
        end
    case 'cart'
        for n=1:tot_k
            kpt_out(n,:)=(PIL_linexpan(kpt_in(n,:),prim_vec'))';
        end
    end
endfunction
