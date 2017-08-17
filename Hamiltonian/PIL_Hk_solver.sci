// **** Purpose ****
// This function is an eigenvalue solver of Hk using parallel run
// **** Variables ****
// [k_point]: nx3, real
// <= all k-points in row vector form
// [lat_vec]: 3x3, real
// <= lattice vectors in row vector form
// [Hr]: tot_k x tot_state x tot_state
// <= the <n(0)|H|m(R)> matrix
// [uc_index]: n x 3, integer
// <= unit cell index of R in terms of n1*a1+n2*a2+n3*a3
// [uc_deg]: n x 1, integer or [], default=[] (optional)
// <= degeneracy of each uc_index, if [], then assume 1 for all. 
// [n_thread]: nx1, integer / [], default=[] (optional)
// <= number of threads for parallel calculation, if 0, use all cores.
//    if [], use serial. default=[] 
// [verb]: 1x1, string, 'on' / 'off'
// <= if verb='on', information of calculation will show up.  
// << Output >> 
// * output variables are optional. If you output less variables, 
//   you can save more memory. 
// [Ek]: tot_state x tot_k, real
// => eigenvalues of each k-points
// [Vk]: tot_state x tot_state*tot_k, complex 
// => eigenvactor matrix of each k-point, in 2D form
// [Hk]: tot_state x tot_state*tot_k, in 2D sparse-form
// => Hamiltonian of each k-point   
// **** Version ****
// 03/23/2016 first built, integrate serial and parallel version
// 03/25/2016 change Hr to 2D form to speed up and make it compatiable 
//            with parallel_run 
// 03/29/2016 add optional output to handle big calculation cases.
// 04/16/2016 set Vk and Hk to native 2D matrix form
// **** Comment ****
// So far, parallel calculaton only supports Scilab 5.x in linux

function [Ek,Vk,Hk]=PIL_Hk_solver(k_point,lat_vec,Hr,uc_index,uc_deg,n_thread,verb)  
    [lhs,rhs]=argn();
    select rhs
    case 4
        uc_deg=[];
        n_thread=[];
        verb='on'
    case 5
        n_thread=[];
        verb='on'    
    end
    // check Hk format
    tot_uc=length(uc_index(:,1));
    if ndims(Hr)==3 then
        tot_state=length(Hr(:,1,1));
        Hr=matrix(Hr,tot_state,tot_state*tot_uc);
    end

    tot_k=length(k_point(:,1));
    tot_state=length(Hr(:,1));
    if verb=='on' then
        mprintf('      PIL_Hk_solver: %d k-points in total\n',tot_k);
    end

    if tot_k >=10 then
        mod_val=ceil(tot_k/10);
    else 
        mod_val=1;
    end

    if n_thread==[] then  // serial 
        if verb=='on' then
            mprintf('        number of thread= serial\n\n');
        end

        select lhs
        case 3
            Ek=zeros(tot_state,tot_k);
            Vk=zeros(tot_state,tot_state*tot_k);
            Hk=spzeros(tot_state,tot_state*tot_k);
            tic(); t1=0;
            for n=1:tot_k
                [hk,err]=PIL_Hk_R(k_point(n,:),lat_vec,uc_index,uc_deg,Hr);
                [V,D]=spec(hk);
                Ek(:,n)=real(diag(D));
                Vk(:,(n-1)*tot_state+1:n*tot_state)=V;
                Hk(:,(n-1)*tot_state+1:n*tot_state)=sparse(hk);
                if verb=='on' then
                    if pmodulo(n,mod_val)==0 then
                        t2=toc();
                        mprintf('        %4d k-points calculated,'..
                        +' time=%f\n',n,t2-t1)
                        t1=t2;
                    end
                end
            end
        case 2
            Ek=zeros(tot_state,tot_k)
            Vk=zeros(tot_state,tot_state*tot_k);
            tic(); t1=0;
            for n=1:tot_k
                [hk,err]=PIL_Hk_R(k_point(n,:),lat_vec,uc_index,uc_deg,Hr);
                [V,D]=spec(hk);
                Ek(:,n)=real(diag(D));
                Vk(:,(n-1)*tot_state+1:n*tot_state)=V;
                if verb=='on' then
                    if pmodulo(n,mod_val)==0 then
                        t2=toc();
                        mprintf('        %4d k-points calculated,'..
                        +' time=%f\n',n,t2-t1)
                        t1=t2;
                    end
                end
            end
        case 1
            Ek=zeros(tot_state,tot_k);
            tic(); t1=0;
            for n=1:tot_k
                [hk,err]=PIL_Hk_R(k_point(n,:),lat_vec,uc_index,uc_deg,Hr);
                [V,D]=spec(hk);
                Ek(:,n)=real(diag(D));
                if verb=='on' then
                    if pmodulo(n,mod_val)==0 then
                        t2=toc();
                        mprintf('        %4d k-points calculated,'..
                        +' time=%f\n',n,t2-t1)
                        t1=t2;
                    end
                end
            end
        end
    else  // parallel
        if verb=='on' then
            mprintf('      PIL_Hk_solver: number of thread=  %d\n\n', n_thread);
        end
        function ek_vk_hk=Hk_solver(k)
            [hk,err]=PIL_Hk_R(k',lat_vec,uc_index,uc_deg,Hr);
            [vk,dk]=spec(hk);
            select lhs
            case 3
                ek_vk_hk=zeros(tot_state+4*tot_state^2,1);
                ek_vk_hk(1:tot_state)=real(diag(dk));
                ek_vk_hk(tot_state+1:tot_state+tot_state^2)=real(vk(:));
                ek_vk_hk(tot_state+tot_state^2+1:tot_state+2*tot_state^2)..
                =imag(vk(:));
                ek_vk_hk(tot_state+2*tot_state^2+1:tot_state+3*tot_state^2)..
                =real(hk(:));
                ek_vk_hk(tot_state+3*tot_state^2+1:tot_state+4*tot_state^2)..
                =imag(hk(:));
            case 2
                ek_vk_hk=zeros(tot_state+2*tot_state^2,1);
                ek_vk_hk(1:tot_state)=real(diag(dk));
                ek_vk_hk(tot_state+1:tot_state+tot_state^2)=real(vk(:));
                ek_vk_hk(tot_state+tot_state^2+1:tot_state+2*tot_state^2)..
                =imag(vk(:));
            case 1
                ek_vk_hk=zeros(tot_state,1);
                ek_vk_hk(1:tot_state)=real(diag(dk));
            end

        endfunction

        // solve eigenvalue problem
        select lhs
        case 3
            Ek_Vk_Hk=zeros(tot_state+4*tot_state^2,tot_k);
            Ek_Vk_Hk=parallel_run(real(k_point'),"Hk_solver",..
            tot_state+4*tot_state^2,init_param('nb_workers', n_thread));
        case 2
            Ek_Vk_Hk=zeros(tot_state+2*tot_state^2,tot_k);
            Ek_Vk_Hk=parallel_run(real(k_point'),"Hk_solver",..
            tot_state+2*tot_state^2,init_param('nb_workers', n_thread));
        case 1
            Ek_Vk_Hk=zeros(tot_state,tot_k);
            Ek_Vk_Hk=parallel_run(real(k_point'),"Hk_solver",..
            tot_state,init_param('nb_workers', n_thread));
        end

        // reshape to normal form
        Ek=Ek_Vk_Hk(1:tot_state,:);
        if lhs==2 | lhs==3
            Vk=..
            matrix(Ek_Vk_Hk(tot_state+1:tot_state+tot_state^2,:),..
            tot_state,tot_state*tot_k)..
            +%i*matrix(..
            Ek_Vk_Hk(tot_state+tot_state^2+1:tot_state+2*tot_state^2,:),..
            tot_state,tot_state*tot_k);
        end
        if lhs==3
            Hk=..
            sparse(matrix(..
            Ek_Vk_Hk(tot_state+2*tot_state^2+1:tot_state+3*tot_state^2,:),..
            tot_state,tot_state*tot_k)..
            +%i*matrix(..
            Ek_Vk_Hk(tot_state+3*tot_state^2+1:tot_state+4*tot_state^2,:),..
            tot_state,tot_state*tot_k));
        end
    end
endfunction

