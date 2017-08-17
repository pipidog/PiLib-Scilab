// **** Purpose ****
// PiLab self-consistent solver (level 2)
// **** variables ****
// ==== << PiLab inputs >> ====
// [scc.HubU]: n x 2, int / empty
// <= U for each state, [state_label, U] or blank 
// [scc.Charge]: 1x total state        
// <= charge of each state, 1x total state
// [scc.Mixing]: 1x1, 0~1 real
// <= mixing parameter, 0~1
// [scc.Iteration]: 1x1, int
// <= maximal number of iterations
// [scc.Converge]: 1x1, < 1e-1 real
// <= convergence criterion, real, at least < 0.1
// [scc.Mesh]: 1x1 / 1x2 / 3x3, int
// <= k-space mesh for calculating Ef, large for metal
// [scc.Temperature]: 1x1, real, default=100
// <= temperature for searching Ef, large for insulator
// ==== << PiLab outputs >> ====
// [scc.E_Fermi]: 1x1, real
// => the chemical potential, may not be accuate if scc.Mesh is small
// [scc.E_gap]: 1x1, real
// => the band gap, may not be accuate if scc.Mesh is small
// [scc.DM_out]: n x 3, real, t-sp
// => the self-consistent density matrix
// [scc.U_mat]: n x 3, real, t-sp
// => the self-consistent Hubbard potential
// **** Version ****
// 05/26/2014 first built
// 05/12/2015 change reload process
// 02/09/2016 remove H_onsite
// 04/16/2016 remove buffer files, use new k-point functions
// **** Comment ****
// 1. This code generates:
//  1).Fermi level (scc.E_Fermi(1))
//  2).Energy Gap (scc.E_gap(1))
//  3).Output density matrix (scc.DM_out(:,:), sparse)
//  4).Hubbard Matrix (scc.U_mat(:,:), sparse)

function PiLab_scc(project_name)
    disp('{scc}: starting calculation ...');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{scc}: loading variables ...');
    PiLab_loader(project_name,'scc','all','trim');
    load(project_name+'_scc.sod'); 
    load(project_name+'_lat.sod');
    load(project_name+'_hop.sod');

    if isfield(scc,'DM_out') then
        if length(scc.DM_out(:,1))~=length(hop.state_info(:,1)) then
            PiLab_loader(project_name,'scc','user','trim');
            load(project_name+'_scc.sod');
            disp('   Note: Previous scc.DM_out automatically erased!')
        else
            disp('   Note: Previous scc.DM_out automatically loaded!')
        end
    end
    if length(scc.HubU)==0 then
        scc.HubU=[1,0];
    end
        
    // check variables =================================================
    disp('{scc}: checking variables ...')
    check_var=(length(scc.HubU(:,1))>=1 & length(scc.HubU(1,:))==2...
    & length(scc.HubU(:,1))<=length(hop.state_info(:,1)));
    if check_var~=%t then
        disp('Error: PiLab_scc, scc.HubU has wrong input!');
        abort;
    end
    check_var=(length(scc.Charge(1,:))==length(hop.state_info(:,1)));
    if check_var~=%t then
        disp('Error: PiLab_scc, scc.Charge has incorrect state number!');
        abort;
    end
    check_var=(scc.Mixing>=0 & scc.Mixing<=1);
    if check_var~=%t then
        disp('Error: PiLab_scc, scc.Mixing must be >=0 & <=1!');
        abort;
    end
    check_var=(scc.Iteration>=0 & fix(scc.Iteration)-scc.Iteration==0);
    if check_var~=%t then
        disp('Error: PiLab_scc, scc.Iteration must >0 & integer!');
        abort;
    end
    check_var=(scc.Converge>=0 & scc.Converge <10^-1);
    if check_var~=%t then
        disp('Error: PiLab_scc, scc.Converge must >0 & <0.1!');
        abort;
    end
    check_var=(length(scc.Mesh)==length(lat.LatVec(1,:)) & sum(scc.Mesh==0)==0 );
    if check_var~=%t then
        disp('Error: PiLab_scc, scc.Mesh has wrong dimension or 0 mesh!');
        abort;
    end
    check_var=(scc.Temperature >=0 );
    if check_var~=%t then
        disp('Error: PiLab_scc, scc.Temperature must >=0!');
        abort;
    end

    // core part ========================================================
    disp('{scc}: running core part ...');
    
    // define necessary parameters
    tot_k=prod(scc.Mesh);
    tot_state=length(scc.Charge);
    tot_e=sum(scc.Charge)*tot_k;

    // generate initial DM_uc_in ---------------------------------------
    disp('  => generating initial density matrix')
    if isfield(scc,'DM_out') & isfield(scc,'U_mat') then
        DM_in=scc.DM_out;
    else
        DM_in=diag(scc.Charge);
    end

    // generate all k-points 
    k_point=PIL_k_mesh(lat.rec_vec,scc.Mesh);
    k_point=PIL_vec_3d(k_point);

    // starting self-consistent calculation
    disp('  => beginning self-consistent calculations')
    
    convergence=100;
    count=0;
    while convergence > scc.Converge
        count=count+1;
        disp('     << iteration '+string(count)+'>>'); 

        // calculate U_mat
        scc.U_mat=zeros(tot_state,tot_state);
        for n=1:length(scc.HubU(:,1))
            scc.U_mat(scc.HubU(n,1),scc.HubU(n,1))...
            =scc.HubU(n,2)*(1/2-DM_in(scc.HubU(n,1),scc.HubU(n,1)));
        end

        // calculate Ef
        // convert hop file to conventional H(R) form 
        [uc_index,H_R]=PIL_H_R(lat.surr_site,hop.state_info,..
        hop.onsite_E+scc.U_mat+hop.LS_mat,hop.hop_mat);
        
        // calculation k-point information
        [Ek,Vk]=PIL_Hk_solver(k_point,lat.LatVec,H_R,uc_index);        
        
        // calculate E_fermi
        scc.E_Fermi=PIL_Ef(Ek(:),tot_e,scc.Temperature);

        // calculate degenerate weighting
        state_weight=tot_e/length(find(Ek(:) <= scc.E_Fermi));

        //integrate all k-space to get new onsite density 
        scc.DM_out=zeros(tot_state,tot_state);        
        for n=1:tot_k
            tot_e_k=length(find(Ek(:,n) <= scc.E_Fermi))*state_weight;
            // read eigenvectors
            scc.DM_out=scc.DM_out...
            +PIL_DM_gen(Ek(:,n),0.01,scc.E_Fermi,..
            tot_e_k,Vk(:,(n-1)*tot_state+1:n*tot_state));
        end
        
        // take average to just one unitcell and erase off-diag terms
        scc.DM_out=real(diag(diag(scc.DM_out)))/tot_k;
        tot_e_diff=abs(sum(diag(scc.DM_out))-sum(scc.Charge))/sum(scc.Charge)
        if tot_e_diff >= 10^-3 then
            disp('Warning: PiLab_scc, total charge difference='...
            +string(tot_e_diff*100)+'% > 0.1%');
        end

        // check convergence
        convergence=max(abs(scc.DM_out-DM_in));
        disp('   Ef='+string(scc.E_Fermi)+' ; state_weight='...
        +string(state_weight)+' ; convergence='+string(convergence));
        DM_in=(1-scc.Mixing)*DM_in+scc.Mixing*scc.DM_out;
    end
    
    // write SCF information
    if convergence > scc.Converge then
        disp('   Exceed Maximal Iterations!');
    else
        disp('   Self-consistency reached!');
    end

    // Band Gap & DOS calculation
    Ek=gsort(Ek(:),'g','i');
    Ef_level=max(find(Ek <= scc.E_Fermi));
    scc.E_gap=Ek(Ef_level+1)-Ek(Ef_level);
    disp('   Band Gap='+string(scc.E_gap));

    // output information ==============================================
    disp('{scc}: output information ...')
    fid(1)=mopen(project_name+'_scc.plb','a+');
    PIL_print_mat('scc.E_Fermi, @f:f, the Fermi level',scc.E_Fermi,'r',fid(1));
    PIL_print_mat('scc.E_gap, @f:f, the band gap',scc.E_gap,'r',fid(1));
    PIL_print_mat('scc.DM_out, @f:ts, the output density matrix'...
    ,sparse(scc.DM_out),'sp',fid(1));
    PIL_print_mat('scc.U_mat, @f:ts, the Hubbard potential matrix'....
    ,sparse(scc.U_mat),'sp',fid(1));
    mclose(fid(1));

    // finishing program ===============================================
    save(project_name+'_scc.sod','scc');
    disp('{scc}: finishing calculation ...');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction 
