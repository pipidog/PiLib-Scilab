// **** Purpose ****
// {kif_esc}: energy surface cut (perform a energy surface cut in 2D BZ)
// **** Variables ****
//    ==== << PiLab inputs >> ====
//    [kif_esc.Command]: 1x1, string
//    <= command to call scilab, e.g:
//    'scilab5 -nwni -f @input_name > @output_name &'
//    [kif_esc.EVal]: 1x1, real
//    <= Energy of the energy surface cut
//    [kif_esc.EWin]: 1x1, real
//    <= Energy window of the energy cut
//    [kif_esc.StateProj]:1xn, integer
//    <= projection of the energy cut on each state.
//    For states selected in that energy window, PiLab can evaluate their
//    weighting on the assigned states. 
//    use negative integer as divider to specify each proejction. 
//    [kif_esc.Thread]: 1x1, integer
//    <= how many thread to run
//    ==== << PiLab outputs >> ====
//    [kif_esc.E_index]: tot_k x tot_proj+1
//    <= [k_point, E_index, n-th weight ]
// **** Version ****
// 05/06/2016: first built
// **** Comment ****

function PiLab_kif_esc(project_name)
    disp('{kif_esc}: starting calculation ...');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{kif_esc}: loading variables ...');
    PiLab_loader(project_name,'kif_esc','user','trim');
    load(project_name+'_kif_esc.sod');
    load(project_name+'_kif.sod');
    disp('  => all variables loaded');

    // variable check ===================================================
    disp('{kif_esc}: checking variables ...');
    // check variables -----------------------------------------------------
    if kif.KptGen~='grid' then
        disp('Error: PiLab_kif_esc, not working for kif.KptGen~=''grid'' case');
        abort
    end
    if length(find(kif.Mesh==1))~=1 then
        disp('Error: PiLab_kif_esc, kif.Mesh is not a 2D mesh!');
        abort
    end
    disp('  => all variables passed !');

    // core Part ========================================================
    disp('{kif_esc} running core part ...');
    work_dir=PIL_dir_path(pwd());
    kif_esc.k_mesh=kif.Mesh(find(kif.Mesh~=1));

    // check kif_esc.StateProj -----------------------------------------
    disp('  => check state projection')
    load(work_dir+'kif/Ek/Ek_1.sod');
    tot_state=length(Ek(:,1));
    if kif_esc.StateProj~=[] then
        if max(kif_esc.StateProj)>tot_state then
            disp('Error: PiLab_kif_esc, unreasonable state projection');
            abort; 
        end
    end

    // create dir ------------------------------------------------------
    disp('  => creating working forlder')
    work_dir=PIL_dir_path(pwd());

    // calculate num of alloc of each thread
    n_alloc=ceil(kif.Allocation/kif_esc.Thread)*ones(1,kif_esc.Thread);
    n_alloc($)=sum(n_alloc(1:$))-sum(n_alloc(1:$-1)); 
    if n_alloc($) <=0  then
        disp('Error: PiLab_kif_esc, negative n_alloc!');
        disp('       Use less thread or make kif.Allocation/kif_esc.Thread=integer!');
        abort;
    end

    // prepare input files ---------------------------------------------
    disp('  => preparing distributed input files')
    mkdir(work_dir+'kif_esc');
    for n=1:kif_esc.Thread
        fid=mopen(work_dir+'kif_esc/input_'+string(n)+'.sce','w');

        mfprintf(fid,'clear;stacksize(''max'');exec(PiLib);\n');
        mfprintf(fid,'load('''+work_dir+project_name+'_kif.sod'');\n');
        mfprintf(fid,'k_start_idx=kif.k_alloc(:,2);\n');
        mfprintf(fid,'clear kif;\n');
        mfprintf(fid,'load('''+work_dir+project_name+'_kif_esc.sod'');\n');
        mfprintf(fid,'E_index=[];\n');
        mfprintf(fid,'for n=%d:%d\n',sum(n_alloc(1:n))-n_alloc(n)+1,sum(n_alloc(1:n)));
        mfprintf(fid,'  load('''+work_dir+'kif/Ek/Ek_''+string(n)+''.sod'');\n');
        mfprintf(fid,'  load('''+work_dir+'kif/Vk/Vk_''+string(n)+''.sod'');\n');
        mfprintf(fid,'  E_index_tmp=PIL_plb_esc_cal(k_start_idx(n),Ek,Vk,kif_esc);\n');
        mfprintf(fid,'  E_index=cat(1,E_index,E_index_tmp);\n');
        mfprintf(fid,'  clear Ek Vk;\n');
        mfprintf(fid,'end\n');    
        mfprintf(fid,'save('''+..
        work_dir+'kif_esc/E_index_'+string(n)+'.sod'',''E_index'');\n');
        mclose(fid);
    end

    // batch submit ----------------------------------------------------
    disp('  => batch submit jobs');
    // reconstruct call command
    // full command is achieved by setting call_command(2) and call_command(4)
    tmp=strsplit(kif.Command);
    fn_index=find(tmp=='@');
    call_command=[strcat(tmp(1:fn_index(1)-1)),'',..
    strcat(tmp(fn_index(1)+11:fn_index(2)-1)),'',..
    strcat(tmp(fn_index(2)+12:$))]
    for n=1:kif_esc.Thread
        // input file
        call_command(2)=work_dir+'kif_esc/input_'+string(n)+'.sce';
        // output file
        call_command(4)=work_dir+'kif_esc/input_'+string(n)+'.log';
        unix(strcat(call_command));
    end

    // detecting job completion ----------------------------------------
    disp('  => detecting job completion');
    E_check=0;
    c2=clock();
    while E_check==0
        E_check=1;
        E_dir=dir(work_dir+'kif_esc/')
        for n=1:kif_esc.Thread
            if grep(E_dir(2),'E_index_'+string(n)+'.sod')==[] then
                E_check=0;
                sleep(20*1e+3);
                break;
            end
        end 
    end
    disp('      all jobs completed ');
    disp('      time elapse= '+string(etime(clock(),c2))+' seconds');

    // collect all data ------------------------------------------------
    disp('  => collecting calculated results')
    kif_esc.E_index=[];
    for n=1:kif_esc.Thread
        load(work_dir+'kif_esc/E_index_'+string(n)+'.sod');
        kif_esc.E_index=cat(1,kif_esc.E_index,E_index);
        clear E_index;
    end

    // combine the same k-points ---------------------------------------
    disp('  => calculating k-weighting');
    tmp=size(kif_esc.E_index)
    kif_esc.k_weight=zeros(tmp(1),tmp(2)-1);
    count=1;
    kif_esc.k_weight(1,:)=kif_esc.E_index(1,[1,3:tmp(2)])
    for n=2:length(kif_esc.E_index(:,1))
        if kif_esc.E_index(n,1)~=kif_esc.E_index(n-1,1) then
            count=count+1;
            kif_esc.k_weight(count,1)=kif_esc.E_index(n,1);
        end
        kif_esc.k_weight(count,2:$)=kif_esc.k_weight(count,2:$)..
        +kif_esc.E_index(n,3:$);
    end
    kif_esc.k_weight=kif_esc.k_weight(find(kif_esc.k_weight(:,1)~=0),:);
    
    // delete temporary folder -----------------------------------------
    disp('  => removing scratch folder')
    removedir(work_dir+'kif_esc');
    
    // output information ==============================================
    disp('{kif_esc} output information ...');
    fid(1)=mopen(project_name+'_kif_esc.plb','a+');
    PIL_print_mat('kif_esc.k_mesh, @f:f, 2D BZ k-mesh', kif_esc.k_mesh,..
    'i',fid(1),'off'); 
    PIL_print_mat('kif_esc.E_index, @f:f, selected eigenvalues,'+..
    ' [k-point, E_index, n-th weight]', kif_esc.E_index,'r',fid(1),'off');
    PIL_print_mat('kif_esc.k_weight, @f:f, [k-point, n-th weight]',..
    kif_esc.k_weight,'r',fid(1),'off');
    mclose(fid(1));

    // finishing program ===============================================
    save(project_name+'_kif_esc.sod','kif_esc');
    disp('{kif_esc}: finishing {kif_esc} calculation ...');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction
