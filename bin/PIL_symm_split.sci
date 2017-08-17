// **** Purpose ****
// This function can generate a small field that breaks crystal field, 
// time-reversal, and parity symmetries of the given "hop" from PiLab.  
// **** Variables ****
// [hop]: structures
// <= the output vairable of hop in PiLab
// [CF_split],[TR_split],[PR_split]: 1*1, real
// <= strength of the split
// [V_pert]: N*N, real
// => perturation potential
// **** Version ****
// 04/09/2015 first built
// **** Comment ****
// 1. Chern and Z2 reuires symmetry breaking to make sure there is not
//    unwanted degeneracy. Note that, for Z2, time-reversal breaking has
//    to be kept! 
// 2. CF_split will generate the following splitting for each orbital 
//    linspace(-CF_split,+CF_split,tot_suborbital)
// 3. TR_split will generate the following splitting for each spin
//    linspace(-TR_split,+TR_split,2)
// 4. PR_split will generate the following splitting for each sublattice
//    linspace(-PR_split,+PR_split,tot_sublattice) 

function V_pert=PIL_symm_split(hop,CF_split,TR_split,PR_split)
    tot_state=2*(2*hop.SiteOrb(:,2)+1);
    tot_site=max(hop.SiteOrb(:,1));
    
    // generate TR & CR spliting table for 32 cubic suboribtals
    orb_split=list()
    
    orb_split(1)=[0] // s-orb
    for n=2:4 // p,d,f-orb
        orb_split(n)=linspace(-CF_split,+CF_split,2*(n-1)+1)
    end    
    site_split=linspace(-PR_split,+PR_split,tot_site);
    
    // generate V_pert in cubic harmonics without state selecting
    V_pert=[];
    L_notation=['s','p','d','f'];
    U_c_out=[];
    for n=1:length(hop.SiteOrb(:,1))
        S=hop.SiteOrb(n,1);
        L=hop.SiteOrb(n,2);
        orb_potential=orb_split(L+1)+site_split(S)*ones(1,2*L+1);
        // spin-dn
        V_pert=cat(2,V_pert,orb_potential-TR_split);
        // spin-up
        V_pert=cat(2,V_pert,orb_potential+TR_split);
        // basis transformation matrix
        [M_out,U_in_out]=PIL_basis_trans(eye(2*(2*L+1),2*(2*L+1))...
        ,L_notation(L+1),'c',hop.Basis);
        U_c_out=PIL_dirsum(U_c_out,U_in_out);
    end
    V_pert=clean(U_c_out*diag(V_pert)*U_c_out');
    V_pert=V_pert(hop.SelState,hop.SelState);
    // check hermitianess of V_pert
    if sum(abs(V_pert-V_pert')) >=10-5 then
        disp('Error: PIL_symm_split, V_pert is not hermitian!');
        abort
    else
       V_pert=(V_pert+V_pert')/2;
   end 
endfunction
  
