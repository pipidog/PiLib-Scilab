// **** Purpose ****
// This function will check whether your Hamiltonian (ham) has TR symmetry
// **** Variables ****
// [ham]: structure
// <= ham variable obtained from ham.
// [k_point]: 1x3, real, optional, default=[3,5,7]
// <= specify a k-point, then the code will check if +k and -k can be
//    linked by TR operator. Optional, defult use: [3,5,7] 
// [criterion: 1x1, real, optiona, default=1e-4
// <= if TR_mat-eye(tot_state,tot_state) has any value larger than 
//    this one, TR is considered non-exist.
// [TR_check]: 1x1, int, 0/1
// => if 0, TR test failed, if 1, TR test passed. 
// [TR_mat]: tot_state x tot_state, real
// => abs(V_TR(k)'*V(-k)), If TR exist, it should be an identity matrix. 
//    (for values small than 1e-6 will be cleaned)
// **** Version ****
// 04/27/2016 first built
// **** Comment ****

function [TR_check,TR_mat]=PIL_TR_check(lat_vec,Hr_mat,uc_index,uc_deg,..
    TR_pair,k_point,criterion)
    [lhs,rhs]=argn()
    select rhs
    case 2        
        k_point=[3,5,7];
        criterion=1e-4
    case 3
        criterion=1e-4
    end
    tot_state=length(Hr_mat(:,1));

    // solve eigen problem
    [Ek_p,Vk_p]=PIL_Hk_solver( k_point,lat_vec,Hr_mat,uc_index,uc_deg,[],'off');
    [Ek_n,Vk_n]=PIL_Hk_solver(-k_point,lat_vec,Hr_mat,uc_index,uc_deg,[],'off');

    // perform TR_check
    [Vk_p_TR]=PIL_TR_op(Vk_p,TR_pair);
    TR_mat=clean(abs(Vk_p_TR'*Vk_n),1e-6);
    if max((TR_mat-eye(tot_state,tot_state)))>=criterion then
        TR_check=0;
    else
        TR_check=1;
    end
endfunction
