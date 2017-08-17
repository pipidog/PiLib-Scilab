// **** Purpose ****
// For a given primitive cell, one can always generates its coventional 
// by assigning three new axes. These axes must be able to presented 
// in terms of the original primitive vectros. 
// **** Variables ****
// [pc_vec]: 3x3, real
// <= the primitive lattice row vectors.
// [pc_sublat]: nx3, real / []
// <= sublattice of primitive cell in cartisian coordinate.   
// [cc_vec]: 3x3, real
// <= the conventional lattice vectors
// [cc_sublat]: nx10, real
// => information of the sublattice in the conventional cell
//    [lattice index in pc(b,n1,n2,n3),x,y,z, expansion coeff in cc_vec]
//    expansion coeff in cc_vec must all range [0,1).   
// **** Version ****
// Jan 28, 2016: first built
// **** Comment ****
// 1. In praticial use, you may just want cc_sublat with xyz position left
//    only, if so, just use cc_sublat=cc_sublat(:,5:7) 

function [cc_sublat]=PIL_conv_cell_vec(pc_vec,pc_sublat,cc_vec)

    // common varables
    tot_pcsub=length(pc_sublat(:,1));
    pc_subred=PIL_red_cart_conv(pc_vec,pc_sublat,'cart')
    scan_order=max(abs(cc_vec))+1;
        
    // constructure cc_vec
    for n=1:3
        cc_vec(n,:)=cc_vec(n,:)*pc_vec
    end
    
    // list all surrounding atoms 
    pc_loop=PIL_nest_loop(repmat([-scan_order,+scan_order],3,1));
    site_list=zeros(length(pc_loop(:,1))*tot_pcsub,7)
    count=0;
    for n=1:length(pc_loop(:,1))
        for m=1:tot_pcsub
            count=count+1;
            site_list(count,1:4)=[m,pc_loop(n,:)];
            site_list(count,5:7)=..
            PIL_lat_index(site_list(count,1:4),pc_vec,pc_sublat);
        end
    end

    // searching for in-supercell atoms
    // cc_sublat=[lattice index in pc, x,y,z, expansion in cc]
    count=0;
    cc_vec_norm=[norm(cc_vec(1,:)),norm(cc_vec(2,:)),norm(cc_vec(3,:))];
    for n=1:length(site_list(:,1))
        exp_coff=(PIL_linexpan(site_list(n,5:7),cc_vec'))';
        if prod(1-exp_coff > 1e-3)==1 & prod(exp_coff >=-1e-3)==1 
            count=count+1;
            cc_sublat(count,:)=[site_list(n,:),exp_coff]
        end
    end
    cc_sublat=gsort(clean(cc_sublat),'lr','i')
endfunction

// example:
//clear; clc; exec(PiLib);
//pc_vec=..
//[   6.300000    0.000000    0.000000
// 4.440000    4.480000    0.000000
//-5.370000   -2.240000    2.430000]
//pc_sublat=..
//[  0.000000    0.000000    0.000000
//3.150000    0.000000    1.215000
//4.478580    1.868160    0.000000
//1.328580    1.868160    1.215000]
//cc_vec=..
// [0 -1 -1
//  -1 0 -1
//  1 1 0]
//[cc_sublat]=PIL_conv_cell_vec(pc_vec,pc_sublat,cc_vec);

