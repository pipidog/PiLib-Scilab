// **** Purpose ****
// This code genertes the flquet H(k) matrix
// **** Variables ****
// [lat],[flq]: structues
// <= the output of lat, scc, flq in PiLab
// [k_point]: 1x3, real
// <= the k_point
// [Hk_flq]: tot_flq_state x tot_flq_state, real
// => the flquet Hamiltonian 
// **** Version ****
// 12/03/2014 first version
// **** Comment ****
// 1. This code is extracted from PiLab_flq_xxx series codes to make 
//    an independent function. 

function Hk_flq=PIL_Hk_flq(k_point,lat,flq,FT_conv)
    [lhs,rhs]=argn();
    if rhs==3 then
        FT_conv='full';
    end

    tot_hop_state=length(flq.state_info(:,1))/(2*flq.Order+1);
    Hk_state_info=cat(2,[1:tot_hop_state]'...
    ,flq.state_info(1:tot_hop_state,3:$));

    // generate Hk for all order
    Hk_sub=list();
    for p=1:flq.Order+1
        Hk_sub(p)=PIL_Hk_gen(k_point,lat.surr_site,Hk_state_info...
        ,flq.H_onsite(p),flq.hop_mat(p),FT_conv,'off');
    end

    // construct flq Hk by linking different order
    tot_flq_state=(2*flq.Order+1)*tot_hop_state;
    Hk_flq=zeros(tot_flq_state,tot_flq_state);
    for p=1:2*flq.Order+1
        for q=p:p+flq.Order
            r_range=(p-1)*tot_hop_state+1:p*tot_hop_state;
            c_range=(q-1)*tot_hop_state+1:q*tot_hop_state;
            if q==p
                Hk_flq(r_range,c_range)=Hk_sub(q-p+1)...
                -((p-1-flq.Order)*flq.Frequency*eye(Hk_sub(q-p+1)));
            elseif q <=2*flq.Order+1
                Hk_flq(r_range,c_range)=Hk_sub(q-p+1);
            else
                break;
            end
        end
    end
    Hk_flq=triu(Hk_flq);
    Hk_flq=Hk_flq+Hk_flq'-diag(diag(Hk_flq));     
endfunction
