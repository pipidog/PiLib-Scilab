// **** Purpose ****
// This file load the PIL_k_solver_dist results into memory 
// **** Variables ****
// [vairable_name]: 3x3, real
// <= description of the variable 
// **** Version ****
// 04/06/2016 first built
// **** Comment ****

function [Ek_dist,Vk_dist,Hk_dist]=PIL_k_dist_loader(work_folder,task_name,n_dist)
    [lhs,rhs]=argn();
    // check work_folder name
    tmp=strsplit(work_folder);
    if tmp($)~='/' | tmp($)~='\' then
        work_folder=work_folder+'/'
    end
    k_folder=work_folder+task_name+'_k/'

    // check dimensions
    load(k_folder+'k_point/k_point_all.sod');
    tot_k=length(k_point_all(:,1));

    load(k_folder+'Ek/Ek_1.sod');
    tot_state=length(Ek(:,1));

    // generate num of k of each allocation
    n_k=ceil(tot_k/n_dist)*ones(1,n_dist);
    if n_dist >1 then
        n_k($)=tot_k-(n_dist-1)*n_k(1);
    else
        n_k=tot_k;
    end

    // read all data
    count=0;
    select lhs
    case 1
        Ek_dist=zeros(tot_state,tot_k);
        for n=1:n_dist
            load(k_folder+'Ek/Ek_'+string(n)+'.sod');
            Ek_dist(:,count+1:count+n_k(n))=Ek
            count=count+n_k(n);
        end 
    case 2
        Ek_dist=zeros(tot_state,tot_k);
        Vk_dist=zeros(tot_state,tot_state,tot_k);
        for n=1:n_dist
            load(k_folder+'Ek/Ek_'+string(n)+'.sod');
            load(k_folder+'Vk/Vk_'+string(n)+'.sod');
            Ek_dist(:,count+1:count+n_k(n))=Ek
            Vk_dist(:,:,count+1:count+n_k(n))=Vk
            count=count+n_k(n);
        end 
    case 3
        Ek_dist=zeros(tot_state,tot_k);
        Vk_dist=zeros(tot_state,tot_state,tot_k);
        Hk_dist=zeros(tot_state,tot_state,tot_k);
        for n=1:n_dist
            load(k_folder+'Ek/Ek_'+string(n)+'.sod');
            load(k_folder+'Vk/Vk_'+string(n)+'.sod');
            load(k_folder+'Hk/Hk_'+string(n)+'.sod');
            Ek_dist(:,count+1:count+n_k(n))=Ek
            Vk_dist(:,:,count+1:count+n_k(n))=Vk
            Hk_dist(:,:,count+1:count+n_k(n))=Hk
            count=count+n_k(n);
        end 
    else
        disp('Error: PIL_k_dist_loader, must assign output variables');
        abort
    end
endfunction
