// **** Purpose ****
// This function generates the hopping matrix of a given crystal 
// structure and
// Slaster-Koster parameters 
// **** Variables ****
// [surr_site]: list(m) -> n x 9, real
// <= this is the output variable of PIL_uc_nb
// [site_orb]: n x 3, integer
// <= the variable specify the orbital of each site.
// [SK_parameter]: n x 7, real
// <= the SK_parameter [orb1, orb2, nn_order, ts, tp, td, tf]
// [hop_fliter]: 1x1, positive real, default=10^-6
// <= the fliter of small matrix elements in H_hop
// [state_info]: n x 5, integer
// => tells how the states are ordered
//    [state_label, site, n(principle), L(angular), sub-orbital]
// [H_hop]: list(tot_sublatt) x tot_suborb x tot_suborb x tot_surr_site, real
// => the hopping matrix of each nearest neighbor order
// **** Version ****
// 05/15/2014 first built
// 05/24/2014 modify output to include r-vector
// 06/03/2014 fix surr_site becomes list problem
// 06/06/2014 full rewirte, each bond has its own H_hop
// **** Comment ****
//  1.site_orb:
//    e.g: [1,1,0; 3,1,2] means sublattice 1 / 2, has an atomic orbital set with 
//          prinipal (or type) quantum number n=1/1 and angular momentum L=0 / 2
//    Note: the principal number n will not have any esseitial function. This is
//          just a number for you to dinguish two orbitals on the same site. 
//          So, when print state info, it will be shown. 
//  2.SK_parameter:
//    e.g: [1,2,ts,tp,td,tf] means orb1 (i.e the orbital defined in site_orb(1,:))
//          and orb2 has SK_parameter: ts,tp,td,tf
//    to prevent from double assign, this function require SK_parameter(1,n) <=
//    SK_parameter(2,n). 
//    Note: The way we label SK_paramater will not use the information of principal 
//          number n in [site_orb]. 
//  3.more infomation:
//    the hop matrix doen't correspond to real Hamiltonian, it just tells you the coupling
//    between two negibhbor sites. The real Hamiltoina has to consider their position vectors 
//    and impose periodic boundary condiition.
//  4.About the H_hop
//    H_hop means tij(ra,rb)\sum{i,j}C_{i}^{+}(ra)C_{j}(rb), so all the H_hop are just 
//    upper or lower triangle. 
function [state_info,H_hop]=PIL_hop_mat(surr_site,site_orb,SK_parameter...
    ,basis_out,hop_order,hop_fliter)
    // make surr_site integers
    [lhs,rhs]=argn();
    select rhs
    case 4
        hop_order=2;        // default set to 2nd order
        hop_fliter=10^-3;   // default
    case 5
        hop_fliter=10^-3;   // default
    end
    // check input
    if length(site_orb(:,1)) > 1 then
        for n=2:length(site_orb(:,1))
            if prod(site_orb(n,:)==site_orb(n-1,:))==1 then
                disp('Error: PIL_hop_mat, two the same orbitals on the same site is not allowed!');
                abort
            end
        end
    end

    for n=1:length(SK_parameter(:,1))
        if round(real(SK_parameter(n,1))) > round(real(SK_parameter(n,2))) then
            disp('Error: PIL_hop_mat, just specify SK_parameter(n,1) smaller than SK_parameter(n,2) terms !');
            abort;
        end
    end    

    // label how user define the orbital before sort
    site_orb=gsort(cat(2,site_orb,[1:length(site_orb(:,1))]'),'lr','i');

    // define necessary parameters
    tot_sublatt=size(surr_site);
    tot_suborb=2*sum(2*site_orb(:,3)+1);

    // sort complex SK_parameter
    SK_sort_tmp=gsort(round(real(SK_parameter(:,1:3))),'lr','i');
    SK_par_tmp=zeros(SK_parameter);
    for n=1:length(SK_sort_tmp(:,1))
        tmp_ind=find(round(real(SK_parameter(:,1)))==SK_sort_tmp(n,1)...
        & round(real(SK_parameter(:,2)))==SK_sort_tmp(n,2) ...
        & round(real(SK_parameter(:,3)))==SK_sort_tmp(n,3))
        SK_par_tmp(n,:)=SK_parameter(tmp_ind,:);
    end
    SK_parameter=SK_par_tmp;


    // generate all SubOrb table, site_suborb --> 
    // [site, n, L, suborb_index, UserOrb_index]
    // the orders of site_suborb is also the state_index
    site_suborb=[];
    for n=1:length(site_orb(:,1))
        //suborb index (defined by PIL_basis_trans), 
        //Orb_index (order of the orbital)
        select site_orb(n,3)
        case 0 // suborb=1~2 ; SK checker=1,1
            site_suborb=cat(1,site_suborb,cat(2,repmat(site_orb(n,:),2,1)...
            ,[1,2]'));
        case 1 // suborb=3~8 ; SK checker=2~4,2~4
            site_suborb=cat(1,site_suborb,cat(2,repmat(site_orb(n,:),6,1)...
            ,[3:8]'));
        case 2 // suborb=9~18 ; SK checker=5~9,5~9
            site_suborb=cat(1,site_suborb,cat(2,repmat(site_orb(n,:),10,1)...
            ,[9:18]'));
        case 3 // subob=19~32 ; SK checker=10~16,10~16
            site_suborb=cat(1,site_suborb,cat(2,repmat(site_orb(n,:),14,1)...
            ,[19:32]'));
        else
            disp('Error: PIL_hop_mat, site_orb(n,3) can only ranges 0 to 3 !');
            abort
        end
    end
    site_suborb(:,[4:5])=site_suborb(:,[5,4]);
    site_suborb=gsort(site_suborb,'lr','i');


    // site_range is also the state_index of the SurOrb of each site
    site_range=zeros(tot_sublatt,2);
    for n=1:tot_sublatt
        tmp=find(site_suborb(:,1)==n);
        site_range(n,1)=min(tmp);
        site_range(n,2)=max(tmp);
    end

    // basis transformation matrix
    U=[];
    for n=1:length(site_orb(:,1))
        select site_orb(n,3)
        case 0
            [M_out,U_in_out]=PIL_basis_trans(eye(2,2),'s','c',basis_out);
        case 1
            [M_out,U_in_out]=PIL_basis_trans(eye(6,6),'p','c',basis_out);
        case 2
            [M_out,U_in_out]=PIL_basis_trans(eye(10,10),'d','c',basis_out);
        case 3
            [M_out,U_in_out]=PIL_basis_trans(eye(14,14),'f','c',basis_out);
        end
        U=PIL_dirsum(U,U_in_out);
    end

    // key part=========================================================
    H_hop=list();
    for n=1:tot_sublatt // run total sublattice
        H_hop(n)=zeros(tot_suborb,tot_suborb,length(surr_site(n)(:,1))-1);
        s1=round(surr_site(n)(1,3));
        for m=2:length(surr_site(n)(:,1)) // run all its neighbor sites
            s2=round(surr_site(n)(m,3));
            r=surr_site(n)(m,7:9)-surr_site(n)(1,7:9);
            order=round(surr_site(n)(m,1));
            if order <= hop_order then
                for p=site_range(s1,1):site_range(s1,2) // run matrix range
                    for q=site_range(s2,1):site_range(s2,2)
                        suborb1=site_suborb(p,4);
                        suborb2=site_suborb(q,4);

                        SK_index=find(round(real(SK_parameter(:,1)))...
                        ==site_suborb(p,5)...
                        & round(real(SK_parameter(:,2)))==site_suborb(q,5)...
                        & round(real(SK_parameter(:,3)))==order);

                        if length(SK_index)==0 then
                            SK_index=find(round(real(SK_parameter(:,1)))...
                            ==site_suborb(q,5)...
                            & round(real(SK_parameter(:,2)))==site_suborb(p,5)...
                            & round(real(SK_parameter(:,3)))==order);
                        end

                        select length(SK_index)
                        case 0 // no assign case
                            H_hop(n)(p,q,m-1)=PIL_SK_int(suborb1,suborb2...
                            ,r,[0,0,0,0]);
                        case 1 // normal case
                            H_hop(n)(p,q,m-1)=PIL_SK_int(suborb1,suborb2...
                            ,r,SK_parameter(SK_index,4:7));
                        else  // conflict assign case
                            disp('Error: PIL_hop_mat, SK_parameter has conflict assignments!');
                            abort;
                        end
                    end
                end
                H_hop(n)(:,:,m-1)=clean(U*H_hop(n)(:,:,m-1)*U');
            else
                // exit calculate higher order hopping. 
                break
            end
        end
    end
    // =================================================================
    state_info=cat(2,[1:tot_suborb]',site_suborb(:,1:4));
endfunction

