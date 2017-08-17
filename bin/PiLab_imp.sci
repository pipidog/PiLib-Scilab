// **** Purpose ****
// This codes can modify the onsite potential of wannier Hr to mimic 
// the effects of simple impurities.
// **** Variables ****
// ==== << PiLab Inputs >> ====
// [imp.HrRead]: 1x1, str, 'wan'/'spl'/'hdr'
// <= specify the Hr to read. 
// [imp.PotState]: 1xn, int
// <= states that onsite fields are applied
// [imp.PotStreng]: 1xn, real
// <= energy shift of each assigned sublattice
// ==== << PiLab Outputs >> ====
// [imp.R0_index]: 1x1, int
// => index of R=0 unit cell in the loaded uc_index
// [imp.pot_orig]: tot_wf x 1, real
// => original onsite potential of each WF.
// [imp.pot_corr]: tot_wf x 1, real
// => onsite potential corrections of each WF. 
// **** Version ****
// 02/24/2016: first built
// **** Comment ****
function PiLab_imp(project_name)
    disp('{imp}: starting calculation ...');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{imp}: loading variables ...');
    PiLab_loader(project_name,'imp','user','trim');
    load(project_name+'_imp.sod');
    load(project_name+'_'+imp.HrRead+'.sod');
    disp('  => all data loaded')

    // check variables  ================================================
    disp('{imp}: checking variables ...')  
    check_var=(imp.HrRead=='wan' | imp.HrRead=='hdr' | imp.HrRead=='spl');
    if check_var~=%t then
        disp('Error: PiLab_imp, imp.HrRead must '...
        +'be ''wan'', ''hdr'', ''spl'' !');
        abort;
    end
    check_var=(length(imp.PotState)==length(imp.PotStreng))
    if check_var~=%t then
        disp('Error: PiLab_imp, imp.PotState and imp.PotStreng '...
        +'must have the same length!');
        abort;
    end
    disp('  => all variables passed')

    // core part =======================================================
    disp('{imp}: running core part ...');
    disp('  => collecting information')
    select imp.HrRead
    case 'wan'
        tot_state=length(wan.state_info(:,1));
        imp.R0_index=PIL_row_find(wan.uc_index(:,1:3),[0,0,0]);
        H0_mat=..
        wan.Hr_mat(:,(imp.R0_index-1)*tot_state+1:imp.R0_index*tot_state)..
        /wan.uc_index(imp.R0_index,4);
    case 'hdr'
        tot_state=length(hdr.state_info(:,1));
        imp.R0_index=PIL_row_find(hdr.uc_index,[0,0,0]);
        H0_mat=..
        hdr.Hr_mat(:,(imp.R0_index-1)*tot_state+1:imp.R0_index*tot_state);
    case 'spl'
        tot_state=length(spl.state_info(:,1));
        imp.R0_index=PIL_row_find(spl.uc_index,[0,0,0])
        H0_mat=..
        spl.Hr_mat(:,(imp.R0_index-1)*tot_state+1:imp.R0_index*tot_state);
    end
    
    H0_mat=full(H0_mat);
    if max(imp.PotState) > tot_state then
        disp('Error: PiLab_imp, imp.PotState '...
        +'has state label that does''t exist !');
        abort;
    end

    disp('  => modifying onsite potentials')
    imp.pot_orig=diag(H0_mat);
    imp.pot_corr=zeros(tot_state,1);
    imp.pot_corr(imp.PotState)=imp.PotStreng';

    // output information ==============================================
    disp('{imp}: output information ...')
    fid=mopen(project_name+'_imp.plb','a+'); 
    PIL_print_mat('imp.R0_index, @f:f, R=0 unit cell in the uc index ',..
    imp.R0_index,'i',fid(1));
    PIL_print_mat('imp.pot_orig, @f:f, original onsite potentials'..
    +' of each WFs',imp.pot_orig,'r',fid(1));
    PIL_print_mat('imp.pot_corr, @f:f, onsite potential corrections'..
    +' of each WFs',imp.pot_corr,'r',fid(1));
    mclose(fid(1))

    // finishing program ===============================================
    save(project_name+'_imp.sod','imp');
    disp('{imp}: finishing calculation ...');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction
