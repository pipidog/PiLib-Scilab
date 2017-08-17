// **** Purpose ****
// Given a primitive unit cell and its wannier Hr couling information. 
// Given another dimension reduced slab. This code tells you sites
// that involved the in <n(0)|H|m(R)> calculations in the new slab cell. 
// **** Variables ****
// [pc_vec]: 3x3, real
// <= primitive lattice row vectors
// [pc_sublat]: nx3, real
// <= pc sublattice in x,y,z coordinate
// [slab_vec]: 3x3, real
// <= Lattice vectors of the slab from PIL_slab_str
// [slab_sublat]: nx6, real
// <= sublattice information from PIL_slab_str
// []: 1x1, integer
// <= finite size vector using the pc lattice vector indexes
// [wan_uc_index]: nx4, integer
// <= wan_uc_index from wan
// 
// [slab_pc_index]: nx4, integer
// => reindex slab sublattice in terms of pc index
// [coup_site]: nx12, real
// => sites that couple to R=0 in super lattice based on the wan_uc_index
//    results. [b,n1,n2,n3,B,N1,N2,N3,x,y,z,proj on fin axis]
//    b,n1,n2,n3: lattice index in primitive cell
//    B,N1,N2,N3: lattice index in super cell
// [coup_slab]: mx3, integer
// => the coupled super cell index
// **** Version ****
// 01/25/2016: first built.
// 02/18/2016: add coup_order variable 
// **** Comment ****
// Note:
// 1. here the wannier coordinate must be identical to hdr coordinete.
// If not, use PIL_rot_points to do it first. 
// 2. The variable pc_sublat should only contain the sublattices that
// really involves the wanniersation.  
// Note: <--   important experience 

function [coup_site,coup_slab,slab_pc_index]=PIL_slab_coup(pc_vec,pc_sublat,..
    slab_vec,slab_sublat,wan_uc_index)

    // common parameters
    tot_slabsub=length(slab_sublat(:,1));
    tot_pcsub=length(pc_sublat(:,1));
    
    // reindex slab sublat using primitive cell 
    sublat_list_tmp=PIL_sublat_list(pc_vec,pc_sublat,1);
    slab_pc_index=zeros(tot_slabsub,4)
    for n=1:tot_slabsub
        slab_pc_index(n,:)..
        =PIL_lat_index(slab_sublat(n,2:4),pc_vec,pc_sublat,..
        [],sublat_list_tmp);
    end

    // determine limits of super lattice coupling size 
    range_tmp=zeros(2,3)
    range_tmp(1,:)=min(wan_uc_index(:,1:3),'r')..
    +min(slab_pc_index(:,2:4),'r');
    range_tmp(2,:)=max(wan_uc_index(:,1:3),'r')..
    +max(slab_pc_index(:,2:4),'r');
    range_list=PIL_nest_loop(range_tmp');

    // coup_site:[index in prim lat, index in slab lat,xyz, proj in fin axis]
    axis_unit=slab_vec(3,:)/norm(slab_vec(3,:))
    r_range=[min(slab_sublat(:,5)),max(slab_sublat(:,5))];
    coup_site=zeros(length(range_list(:,1))*tot_pcsub,12);
    count=0
    sublat_list_tmp=PIL_sublat_list(slab_vec,slab_sublat(:,2:4),1);
    tot_range_list=length(range_list(:,1));
    printf('\n');
    for n=1:tot_range_list
        for m=1:tot_pcsub
            r=PIL_lat_index([m,range_list(n,:)],pc_vec,pc_sublat)
            r_dot_u3=r*axis_unit'
            if (r_dot_u3>=r_range(1)-0.2) & (r_dot_u3<=r_range(2)+0.2)
                count=count+1;
                // index in prim_vec
                coup_site(count,1:4)=[m,range_list(n,:)] 
                // cartisan position
                coup_site(count,9:11)=r; 
                // index in slab_vec
                coup_site(count,5:8)=PIL_lat_index(coup_site(count,9:11)..
                ,slab_vec,slab_sublat(:,2:4),[],sublat_list_tmp) 
                // proj in fin axis
                coup_site(count,12)=r_dot_u3 
            end 
        end
        if pmodulo(n,fix(tot_range_list/10))==0 then
            printf('     %6.2f%% processed\n',100*n/tot_range_list);
        end
    end
    coup_site=coup_site(1:count,:)
   
    // eliminate repeated coupled slab cell
    tmp=gsort(coup_site(:,6:8),'lr','i')
    count=0;
    coup_slab=100*ones(tmp);
    coup_slab(1,:)=tmp(1,:);
    for n=2:length(tmp(:,1))
        // check if the same
        if sum(abs(tmp(n,:)-tmp(n-1,:)))~=0 then
            coup_slab(n,:)=tmp(n,:)
        end        
    end
    coup_slab=coup_slab(find(coup_slab(:,1)~=100),:)
endfunction
