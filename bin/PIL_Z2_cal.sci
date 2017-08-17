// **** Purpose ****
// This code calculates the Z2 invariant of a 2D plane of a 3D object
// For 2D object, just input appropriate parameters.
// This code can calculate both Z2 or Floquet Z2 
// **** Variables ****
// [lat],[hop],[scc],[flq]: structues
// <= the output of lat, hop, scc, flq in PiLab
//    if Z2, let [flq]=[], if FZ2, let [hop]=[] & [scc]=[] 
// [b_fix]: 1x2, int
// <= Only read when 3D, i.e lat.Primitive is 3x3
//    b_fix(1): The fixed axis 
//    b_fix(2): The fixed value, 0 for 0 , 1 for %pi
// [b_mesh]: 1x2, int
// <= mesh of the half BZ. eg.[5,5] means you mesh the whole BZ from
//    [-%pi,%pi] into 11x11. (always odd to include TRIM)
// [occ_band}: 1x1, int
// <= how many occupied bands below Ef
// [Z2_val]: 1x1, int
// => The Z2 invariant, can be 0 or 1
// [n_field]: tot_mesh points x 3
// => The n_field, [j1,j2,n_field]
//    n_field=F_field+A1_diff-A2_diff, must be 0 or 1
// **** Version ****
// 04/30/2014 first version
// 06/06/2015 combine PIL_Z2_cal & PIL_Z2_flq
// **** Comment ****
// 1. These formulas can be found in JPSJ 76 053702. Note that all of their
// formulas are based on full-periodic Bloch functions (FPBFs) rather
// than cell-periodic Bloch functions (CPBFs). In Comp Phys Comm 183
// 1849, they uses CPBFs, so periodic gauge has to be properly handled.
// 2. This function determins Z2 or FZ2 by checking whether hop, scc, flq
// are empty objects. 

function [Z2_val,n_field]=PIL_Z2_cal(lat,hop,scc,flq,b_fix,b_mesh,occ_band)
    // check task type =================================================
    if flq==[] & hop~=[] & scc~=[] then // hop
        task=1 
    elseif flq~=[] & hop==[] & scc==[]  // flq
        task=2 
    else 
        disp('Error: PIL_Z2_cal, task cannot be specified!');
    end
    // check dimension =================================================
    lat_size=size(lat.Primitive);
    if lat_size(1)==2 then
        lat.recip_vec=PIL_vec_3d(lat.recip_vec);
        lat.recip_vec(3,:)=[0,0,0];
        b_fix=[3,0];
    end

    // define & check paremeters ========================================
    // reorder lat.recip_vec 
    if b_fix(2)~=0 & b_fix(2)~=1 then
        disp('Error:PIL_Z2_cal, b_fix(2) must be 0 or 1');
        abort
    end
    b_vec=lat.recip_vec(find([1:3]~=b_fix(1)),:);

    // tot_state & TR_pair
    tot_k=prod(2*b_mesh+1);
    select task
    case 1
        tot_state=length(hop.state_info(:,1));
        TR_pair=PIL_TR_pair(hop.state_info);
    case 2
        if prod(flq.Phase==[0,0])~=1 & prod(flq.Phase==[0,0,0])~=1
            disp('Error:PIL_Z2_cal, flq.Phase must be all zero!');
            abort;
        end
        tot_state=length(flq.state_info(:,1));
        TR_pair=PIL_TR_pair(flq.state_info);
    end

    // classify k-points ===============================================
    // [kx,ky,kz,j1,j2,zone,k_TR, k+1, k+2]
    // zone: -1 (lower), 0 (TRIM), +1 (upper)
    k_point=zeros(tot_k,9);
    count=0
    for n1=1:2*b_mesh(1)+1
        for n2=1:2*b_mesh(2)+1
            count=count+1;
            j1=((n1-1)/(2*b_mesh(1)))-(1/2);
            j2=((n2-1)/(2*b_mesh(2)))-(1/2);
            k_point(count,1:3)=j1*b_vec(1,:)+j2*b_vec(2,:)...
            +b_fix(2)*lat.recip_vec(b_fix(1),:)/2;
            k_point(count,4:5)=[j1,j2];

            // zone classification
            if j2==-1/2  then // -%pi boundary
                if j1 == -1/2 | j1 == 0 | j1 == +1/2  
                    k_point(count,6)=0;
                elseif j1 < 0 
                    k_point(count,6)=+1;
                elseif j1 > 0
                    k_point(count,6)=-1;
                end
            elseif j2==0 // 0 boundary
                if j1 == -1/2 |  j1 == 0 | j1 == +1/2 
                    k_point(count,6)=0;
                elseif j1 < 0
                    k_point(count,6)=-1;
                elseif j1 > 0
                    k_point(count,6)=+1;
                end
            elseif j2==1/2 // +%pi boundary
                if j1 == -1/2 |  j1 == 0 | j1 == +1/2 
                    k_point(count,6)=0;
                elseif j1 < 0
                    k_point(count,6)=+1;
                elseif j1 > 0
                    k_point(count,6)=-1;
                end
            elseif j2 < 0 // lower plane
                k_point(count,6)=-1;
            elseif j2 > 0 // upper plane
                k_point(count,6)=+1;
            end

            // TR dual
            k_point(count,7)=tot_k-count+1;    

            // generate k+1 & k+2
            if k_point(count,4)~=1/2
                k_point(count,8)=count+2*b_mesh(1)+1;
            end
            if k_point(count,5)~=1/2
                k_point(count,9)=count+1;
            end       
        end
    end

    // impose gauge ====================================================
    // generate eigenstates of B- and B+ -------------------------------
    //disp('PIL_Z2_cal: computing wave functions and guage fixing ...');
    k_band=zeros(tot_state,tot_k);
    k_vec=zeros(tot_state,tot_state,tot_k);

    // B- & B+ not on edge
    B_zone=find(k_point(:,6)==-1 ...
    & (k_point(:,4)~=1/2 & k_point(:,5)~=1/2));
    for n=1:length(B_zone)
        k_ind=B_zone(n);
        k_TR=k_point(k_ind,7);
        select task
        case 1
            Hk=PIL_Hk_gen(k_point(k_ind,1:3),lat.surr_site,hop.state_info...
            ,scc.H_onsite,hop.hop_mat,'full');
        case 2
            Hk=PIL_Hk_flq(k_point(k_ind,1:3),lat,flq,'full');
        end
        [V,D]=spec(Hk);
        // B- area
        k_band(:,k_ind)=diag(D);
        k_vec(:,:,k_ind)=V;
        // check other degeneracy
//        if (min(abs(k_band(1:2:$,k_ind)-k_band(2:2:$,k_ind)))<=1e-7)... 
//            | (min(abs(k_band(3:2:$-1,k_ind)-k_band(2:2:$-1,k_ind)))<=1e-7) then
//            disp('Warning: PIL_Z2_cal, degeneracy beyond TR found at '...
//            +'k_coff=['+string(k_point(k_ind,4))+','+string(k_point(k_ind,5))+']');
//        end
        // B+ area
        k_band(:,k_TR)=k_band(:,k_ind);
        k_vec(:,:,k_TR)=PIL_TR_op(V,TR_pair);
    end

    // B- on edge
    B_zone=find(k_point(:,6)==-1 ...
    & (k_point(:,4)==1/2 | k_point(:,5)==1/2));
    for n=1:length(B_zone)
        k_ind=B_zone(n);
        k_TR=k_point(k_ind,7);
        k_pbc=0;
        if k_point(k_ind,4)==1/2 then
            k_pbc=PIL_row_find(k_point(:,4:5),[-1/2,k_point(k_ind,5)]);
        elseif k_point(k_ind,5)==1/2
            k_pbc=PIL_row_find(k_point(:,4:5),[k_point(k_ind,4),-1/2]);
        end
        // B- area
        k_band(:,k_ind)=k_band(:,k_pbc);
        k_vec(:,:,k_ind)=k_vec(:,:,k_pbc);
        //B+ area
        k_band(:,k_TR)=k_band(:,k_ind);
        k_vec(:,:,k_TR)=PIL_TR_op(k_vec(:,:,k_ind),TR_pair);
    end


    // B0 not on edge
    B_zone=find(k_point(:,6)==0 ...
    & (k_point(:,4)~=1/2 & k_point(:,5)~=1/2));
    odd_state=[1:2:tot_state];
    even_state=[2:2:tot_state]; 
    for n=1:length(B_zone)
        k_ind=B_zone(n);
        select task
        case 1
            Hk=PIL_Hk_gen(k_point(k_ind,1:3),lat.surr_site,hop.state_info...
            ,scc.H_onsite,hop.hop_mat,'full');
        case 2
            Hk=PIL_Hk_flq(k_point(k_ind,1:3),lat,flq,'full');
        end
        [V,D]=spec(Hk);
        k_band(:,k_ind)=diag(D);
        // check degeneracy
//        if min(abs(k_band(1:2:$-2,k_ind)-k_band(3:2:$,k_ind)))<=1e-7 then
//            disp('Warning: PIL_Z2_cal, degeneray beyond TR found at '...
//            +'k_coff=['+string(k_point(k_ind,4))+','+string(k_point(k_ind,5))+']');
//        end
        // odd state
        k_vec(:,odd_state,k_ind)=V(:,odd_state);
        // even state
        V=PIL_TR_op(V,TR_pair);
        k_vec(:,even_state,k_ind)=V(:,odd_state);
    end

    // B0 on edge 
    B_zone=find(k_point(:,6)==0 ...
    & (k_point(:,4)==1/2 | k_point(:,5)==1/2));
    for n=1:length(B_zone)
        k_ind=B_zone(n);
        if k_point(k_ind,4)~=0 & k_point(k_ind,5)~=0 
            k_pbc=PIL_row_find(k_point(:,4:5),[-1/2,-1/2]);
        elseif k_point(k_ind,4)==0
            k_pbc=PIL_row_find(k_point(:,4:5),[0,-k_point(k_ind,5)]);
        elseif k_point(k_ind,5)==0
            k_pbc=PIL_row_find(k_point(:,4:5),[-k_point(k_ind,4),0]);
        else
            disp('Error: PIL_Z2_cal, B0 k_pbc not found!');
            abort
        end
        k_band(:,k_ind)=k_band(:,k_pbc);
        k_vec(:,:,k_ind)=k_vec(:,:,k_pbc);
    end
    // check empty k_vec
    for n=1:tot_k
        if PIL_equal(sum((abs(k_vec(:,:,k_ind))).^2),tot_state)~=%t then
            disp('Error: PIL_Z2_cal, k_vec is not conserved!');
            abort
        end
    end

    // define Berry phases functions ===================================
    function U_val=U_link(basis_ind,k_ind)
        U_val=det(k_vec(:,1:occ_band,k_ind)'...
        *k_vec(:,1:occ_band,k_point(k_ind,7+basis_ind)));
        U_val=inv(abs(U_val))*U_val;
    endfunction

    function F_val=F_field(k_ind)
        k1_ind=k_point(k_ind,8);
        k2_ind=k_point(k_ind,9);
        F_val=log(U_link(1,k_ind)*U_link(2,k1_ind)...
        *inv(U_link(1,k2_ind))*inv(U_link(2,k_ind)));
        if abs(real(F_val)) >= 1e-5 then
            disp('Error: PIL_Z2_cal, F_val is not pure imaginary!');
            abort;
        else
            F_val=imag(F_val)
            // move to main branch
            F_val=F_val-(2*%pi)*round(F_val/(2*%pi));
            if abs(F_val+%pi) < 1e-6 then
                F_val=%pi;
            end
        end
    endfunction

    function A_val=A_field(basis_ind,k_ind)
        A_val=log(U_link(basis_ind,k_ind));
        if abs(real(A_val)) >= 1e-5 then
            disp('Error: PIL_Z2_cal, A_val is not pure imaginary!');
            abort;
        else
            A_val=imag(A_val)
            // move to main branch
            A_val=A_val-(2*%pi)*round(A_val/(2*%pi));
            if abs(A_val+%pi) < 1e-6 then
                A_val=%pi;
            end 
        end
    endfunction

    function A_fd_val=A_fd(basis_fd,basis_ind,k_ind)
        k_fd_ind=k_point(k_ind,7+basis_fd);
        A_fd_val=A_field(basis_ind,k_fd_ind)-A_field(basis_ind,k_ind);
    endfunction
    // calculate n-field ===============================================
    //disp('PIL_Z2_cal: computing n-field ...');
    BZ_zone=find(k_point(:,4)<1/2 & k_point(:,5)<1/2);
    tot_BZ=length(BZ_zone);
    n_field=zeros(tot_BZ,3);
    for n=1:tot_BZ
        k_ind=BZ_zone(n);
        n_field(n,1:2)=k_point(k_ind,4:5);
        n_field(n,3)=(F_field(k_ind)-A_fd(1,2,k_ind)+A_fd(2,1,k_ind))/(2*%pi);
    end
    
    if PIL_equal(n_field(:,3),round(n_field(:,3))) then
        n_field(:,3)=round(n_field(:,3));
    else
        disp('Error:PIL_Z2_cal, n_field are not integers!');
        abort;
    end
    if abs(sum(n_field(:,3))) >= 1e-5  then
        disp('Error:PIL_Z2_cal, sum over all n_field is not zero!');
        abort;
    end
    Z2_val=pmodulo(sum(n_field(find(n_field(:,2)<0),3)),2);
endfunction
