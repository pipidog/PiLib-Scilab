// **** Purpose ****
// finds the surrounding sites of a unitcell up to Nth order
// **** Variables ****
// [primitive]: 3x3 / 2x2 / 1x1, real
// <= the primitive vectros in row
// [sublatt]: nx3 / nx2 / nx1, real
// <= the (x,y,z) cartisian coordinate of sublattices in the unitcell
// [vec_order]: 1x1, integer 
// <= c1*a1+c2*a2+c3*a3, if vec_order=2, then it lists all neighbors
//    within -2 <= [c1,c2,c3] <= +2, must larger than 1
// [NN_criterion]: 1x1, real
// <= criterion of labeling the NN order. if the difference between
//    is small than this valus, they are consider as the same NN order. 
// [NN_order]:1x1, int
// <= neighbor order higer than this will be ignored in surr_site.
// [surr_site]: list(total_sublatt) -> total_nb x 9 x  , real
// => the surrounding sites of that sublattice up to N-th vec_order
//    [nn_order,distant, sublattice label, n1, n2, n3, x, y, z]  
//    r=n1*a1+n2*a2+n3*a3+(sublatt(n,:))=[x,y,z]
// **** Version ****
// 05/01/2014 first built
// 05/20/2014 full rewrite, performance improved, accept low dimension
//            ba_ratio, ca_ratio inputs removed.
// 05/29/2014 fix bug, distant equal set to 10^(-6);
// 06/03/2014 fix bug, surr_site may have different for each sublatt, 
//            so surr_site has been modified to list-type!
// 01/07/2016 add vec_order and NN_criterion. So the degree of search
//            is defined by user. Performance and readability have
//            also much improved.
// **** Comment ****
// 1. This function accepts low dimension input. e.g , if your premitive
//    cell and sublatt are two dimenstion vectors, the code will 
//    generate surr_site table with n3 and z as zero. Use it carefully. 
// 2. One should tune vec_order to make sure the nn_order is correct
//    for high order NNs.
// 3. NN: n-th neighbors
// 4. Don't put too regious value to NN_criterion, it should be around 
//    0.1 to get reasonable results

function [surr_site]=PIL_uc_nb(primitive,sublatt,vec_order,NN_criterion,NN_order);
    // variable check 
    if length(primitive(1,:))~=length(sublatt(1,:)) then
        disp('Error: PIL_uc_nb, dimeisnion inconsistent!');
        abort; 
    end

    if  vec_order <1 then
        disp('Error: PIL_uc_nb, vec_order must be greater than 1')
        abort
    end

    // generate unit cell index 
    dim=length(primitive(1,:));
    tot_sublat=length(sublatt(:,1));
    select dim
    case 3
        loop_index=PIL_nest_loop([-vec_order,vec_order;-vec_order,vec_order;-vec_order,vec_order])
    case 2
        loop_index=PIL_nest_loop([-vec_order,vec_order;-vec_order,vec_order])
    case 1
        loop_index=PIL_nest_loop([-vec_order,vec_order])
    end
    
    // search all sites reside in the super cell
    count=0;
    site_pos=zeros(tot_sublat*(2*vec_order+1)^dim,7);
    for n=1:length(loop_index(:,1))
        rc=loop_index(n,:)*primitive;
        for m=1:tot_sublat
            count=count+1;
            site_pos(count,1:4)=[m,loop_index(n,:)];
            site_pos(count,5:7)=PIL_vec_3d(rc+sublatt(m,:));
        end
    end
    
    // construct surr_site
    tot_site=length(site_pos(:,1));
    surr_site=list();
    for n=1:tot_sublat
        surr_site(n)=zeros(tot_site,9)
        //[nn_order,distant, sublattice label, n1, n2, n3, x, y, z]
        for m=1:tot_site
            d=norm(site_pos(m,5:7)-sublatt(n,:));
            surr_site(n)(m,:)=[0,d,site_pos(m,:)]
        end
        surr_site(n)=PIL_lsort(surr_site(n),'c',[2:9,1],'i')
        
        // label order
        order=0;
        for m=2:tot_site
            if ((surr_site(n)(m,2)-surr_site(n)(m-1,2))) >= NN_criterion then
                order=order+1;
            end
            if order > NN_order then
                // quit for higher order neighbors
                surr_site(n)=surr_site(n)(1:m-1,:);
                break;
            else
                // label neighbor order
                surr_site(n)(m,1)=order;
            end
        end
        
        // reorder based based on sublatt label
        surr_site(n)=PIL_lsort(surr_site(n),'c',[1,3:9,2],'i')
    end
endfunction

// examples of this function: (TaAs structure)
//primitive=..  // data from ab initio
//[  6.305100000    0.000000000    0.000000000;..
//   4.439200000    4.477400000    0.000000000;..
//  -5.372100000   -2.238700000    2.425300000 ]
//sublatt=..
//[0.00000   0.00000   0.00000;..
// 3.15258   0.00000   1.21265;..
// 4.48037   1.86708   0.00000;..
// 1.32785   1.86708   1.21265]
// 
//vec_order=1
//NN_criterion=0.1
