// **** Purpose ****
//  This function calculates the reciprocal lattice vectors of slab structure 
// **** Variables ****
// [lat_vec]: real, 3x3, 
// <= the row lattice vectors where a1 and a2 defines the slab plane. 
// [z_axis]: int, 1/2/3
// <= which axis will be the finite size axis
// [G0]: real, 3x3, 
// => the original reciprocal lattice
// [Gslab]: real, 2x3 or [] if not found
// => the slab reciprocal lattice vectors
// **** Version ****
// 10/18/2016: 1st version
// **** Comment ****

function [Gslab,G0]=PIL_red_BZ_vec(lat_vec,z_axis)
    [lhs,rhs]=argn();
    if rhs==1 then
        z_axis=[];
    end
    G0=PIL_recip_vec(lat_vec);
    select z_axis
    case 2
        Gslab=PIL_recip_vec(diag([1,1e+6,1])*lat_vec);
    case 1
        Gslab=PIL_recip_vec(diag([1e+6,1,1])*lat_vec);
    else
        Gslab=PIL_recip_vec(diag([1,1,1e+6])*lat_vec);
    end

    // check if there is any axis becomes zero
    Gslab_norm=zeros(1,3);
    for n=1:3
        Gslab_norm(n)=norm(Gslab(n,:));
    end
    Gslab=Gslab(find(Gslab_norm>=1e-2),:);
    if length(Gslab(:,1))~=2 then
        disp('Warning: slab basis is not defined');    
    end
endfunction
