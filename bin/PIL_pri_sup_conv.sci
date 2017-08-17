// **** Purpose ****
// Give a primitive cell and a super cell. If their xyz are defined in 
// two different coordinates. This code can convert the primitive cell 
// to the super cell coordinate by specifying three identical points
// in both lattice.
// **** Variables ****
// [eq_site1]: 3x4, integer
// <= the three equal sites in primitive cell in [b,n1,n2,n3] index
// [prim_vec]: 3x3, real
// <= primitive row vectors
// [prim_sub]: nx3, real
// <= primitive sublattice row vectors  
// [eq_site2]: 3x4, integer
// <= the three equal sites in super cell in [b,n1,n2,n3] index
// [sup_vec]: 3x3, real
// <= super row vectors
// [sup_sub]: nx3, real
// <= super sublattice row vectors 
// [prim_vec_out]: 3x3, real
// => the primitive vectors in super cell's xyz coordinate
// [prim_sub_out]: nx3, real
// => the primitive sublattices in super cell's xyz coordinate
// **** Version ****
// Jan 25, 2016: first built
// **** Comment ****
function [prim_vec_out,prim_sub_out]=PIL_pri_sup_conv(eq_site1,prim_vec,..
    prim_sub,eq_site2,sup_vec,sup_sub)

    pt1=zeros(3,3);
    pt2=zeros(3,3);
    for n=1:3
        pt1(n,:)=PIL_lat_index(eq_site1(n,:),prim_vec,prim_sub);
        pt2(n,:)=PIL_lat_index(eq_site2(n,:),sup_vec,sup_sub);
    end
    [M,V_shift,str_diff]=PIL_rot_points(pt1,pt2);

    for n=1:3
        prim_vec_out(n,:)=(M*prim_vec(n,:)')'
    end
    for n=1:length(prim_sub(:,1))
        prim_sub_out(n,:)=(M*prim_sub(n,:)')'
    end
endfunction
