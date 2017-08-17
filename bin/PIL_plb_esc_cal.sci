// **** Purpose ****
// This code is to help calculate PiLab_kif_esc in parallel way. 
// **** Variables ****
// [Ek],[Vk]: data readed from /kif/Ek/Ek_x.sod, /kif/Vk/Vk_x.sod.
// [EVal],[EWin],[StateProj]: kif_esc input data
// [k_idx_start]: the start index of Ek and vk
// **** Version ****
// 05/18/2016 first built
// **** Comment ****


function E_index=PIL_plb_esc_cal(k_idx_start,Ek,Vk,kif_esc)
    EVal=kif_esc.EVal;
    EWin=kif_esc.EWin;
    StateProj=kif_esc.StateProj;
    tot_state=length(Ek(:,1));
    // search target eigenstate and k index 
    [s_ind,k_ind]=find(abs(Ek-EVal)<=EWin); 

    // calculate wieght on projected states
    E_index=[];
    proj_ind=find(StateProj<0);
    if s_ind~=[] then
        if proj_ind==[] then
            w_val=ones(length(s_ind),1);
        else
            w_val=zeros(length(s_ind),length(proj_ind)+1);
            w_val(:,1)=ones(length(s_ind),1);            
            for m=1:length(s_ind)        
                for p=1:length(proj_ind)
                    if p==length(proj_ind) then
                        proj_state=StateProj(proj_ind(p)+1:$);
                    else
                        proj_state=..
                        StateProj(proj_ind(p)+1:proj_ind(p+1)-1);
                    end            
                    w_val(m,p+1)=..
                    sum((abs(Vk(proj_state,..
                    (k_ind(m)-1)*tot_state+s_ind(m)))).^2);
                end
            end
        end
        k_ind=k_ind+(k_idx_start-1)
        E_index=cat(1,E_index,[k_ind',s_ind',w_val]);
    end
    E_index=gsort(E_index,'lr','i');
endfunction
