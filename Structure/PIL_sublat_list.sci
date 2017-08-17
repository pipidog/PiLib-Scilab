// **** Purpose ****
// This function lists all the sublattices within the assigned unit cell
// range.  
// **** Variables ****
// [lat_vec]: 3x3, real
// <= lattice row vectors
// [sublat]: nx3, real
// <= cartisian row vectors of each sublattice
// [vec_order]: 1x1, integer
// <= the unit cell range: n1*a1+n2*a2+n3*a3. If vec_order=1, then 
//    all ni's=-1 ~ +1
// [sublat_list]: tot_uc*tot_sublat x 7
// => sublattices within the assigned ranges. 
//    [uc_index, x,y,z] 
// **** Version ****
// 02/25/2016 first built
// **** Comment ****

function sublat_list=PIL_sublat_list(lat_vec, sublat, vec_order)
    uc_list=PIL_nest_loop([-vec_order,vec_order;..
    -vec_order,vec_order;-vec_order,vec_order])..

    tot_uc=length(uc_list(:,1));
    tot_sublat=length(sublat(:,1));
    sublat_list=zeros(tot_uc*tot_sublat,7)
    for n=1:tot_uc
        sublat_list((n-1)*tot_sublat+1:n*tot_sublat,:)=..
        cat(2,cat(2,[1:tot_sublat]',repmat(uc_list(n,:),tot_sublat,1)),..
        repmat(uc_list(n,:)*lat_vec,tot_sublat,1)+sublat)
    end
endfunction
