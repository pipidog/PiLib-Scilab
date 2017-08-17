// **** Purpose ****
// {kif_sbn}: energy surface cut (perform a energy surface cut in 2D BZ)
// **** Variables ****
//    ==== << PiLab inputs >> ====
//    [kif_sbn.SelBan]: 1xn, integer
//    <= band index of the selected bands
//    [kif_sbn.Verbosity]: 1x1, string, 'more' / 'less'
//    <= whether to print out all eigenvalues. 
//    [kif_sbn.Thread]: 1x1, integer
//    <= how many threads to run
//    ==== << PiLab outputs >> ====
//    [kif_sbn.k_mesh]: 1x3, integer
//    <= 2D kmesh
//    [kif_sbn.k_band]: tot_SelBan * tot_k
//    <= band structure of each band structure
// **** Version ****
// 05/19/2016: first built
// **** Comment ****

function PiLab_kif_sbn(project_name)
    disp('{kif_sbn}: starting calculation ...');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{kif_sbn}: loading variables ...');
    PiLab_loader(project_name,'kif_sbn','user','trim');
    load(project_name+'_kif_sbn.sod');
    load(project_name+'_kif.sod');
    disp('  => all variables loaded');

    // variable check ===================================================
    disp('{kif_sbn}: checking variables ...');
    // check variables -----------------------------------------------------
    if kif.KptGen~='grid' then
        disp('Error: PiLab_kif_sbn, not working for kif.KptGen~=''grid'' case');
        abort
    end
    if length(find(kif.Mesh==1))~=1 then
        disp('Error: PiLab_kif_sbn, kif.Mesh is not a 2D mesh!');
        abort
    end
    disp('  => all variables passed !');

    // core Part ========================================================
    disp('{kif_sbn} running core part ...');
    work_dir=PIL_dir_path(pwd());
    kif_sbn.k_mesh=kif.Mesh(find(kif.Mesh~=1));

    // create dir ------------------------------------------------------
    disp('  => creating working forlder')
    work_dir=PIL_dir_path(pwd());

    // calculate num of alloc of each thread
    n_alloc=ceil(kif.Allocation/kif_sbn.Thread)*ones(1,kif_sbn.Thread);
    n_alloc($)=sum(n_alloc(1:$))-sum(n_alloc(1:$-1)); 
    if n_alloc($) <=0  then
        disp('Error: PiLab_kif_sbn, negative n_alloc!');
        disp('       Use less thread or make kif.Allocation/kif_sbn.Thread=integer!');
        abort;
    end

    // prepare input files ---------------------------------------------
    disp('  => preparing distributed input files')
    mkdir(work_dir+'kif_sbn');
    for n=1:kif_sbn.Thread
        fid=mopen(work_dir+'kif_sbn/input_'+string(n)+'.sce','w');

        mfprintf(fid,'clear; stacksize(''max''); exec(PiLib);\n');
        mfprintf(fid,'load('''+work_dir+project_name+'_kif.sod'');\n');
        mfprintf(fid,'k_start_idx=kif.k_alloc(:,2);\n');
        mfprintf(fid,'clear kif;\n');
        mfprintf(fid,'load('''+work_dir+project_name+'_kif_sbn.sod'');\n');
        mfprintf(fid,'k_band=[];\n');
        mfprintf(fid,'for n=%d:%d\n',sum(n_alloc(1:n))-n_alloc(n)+1,sum(n_alloc(1:n)));
        mfprintf(fid,'  load('''+work_dir+'kif/Ek/Ek_''+string(n)+''.sod'');\n');
        mfprintf(fid,'  k_band=cat(2,k_band,Ek(kif_sbn.SelBan,:))\n');
        mfprintf(fid,'  clear Ek;\n');
        mfprintf(fid,'end\n');    
        mfprintf(fid,'save('''+..
        work_dir+'kif_sbn/k_band_'+string(n)+'.sod'',''k_band'');\n');
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
    for n=1:kif_sbn.Thread
        // input file
        call_command(2)=work_dir+'kif_sbn/input_'+string(n)+'.sce';
        // output file
        call_command(4)=work_dir+'kif_sbn/input_'+string(n)+'.log';
        unix(strcat(call_command));
    end

    // detecting job completion ----------------------------------------
    disp('  => detecting job completion');
    E_check=0;
    c2=clock();
    while E_check==0
        E_check=1;
        E_dir=dir(work_dir+'kif_sbn/')
        for n=1:kif_sbn.Thread
            if grep(E_dir(2),'k_band_'+string(n)+'.sod')==[] then
                E_check=0;
                sleep(5*1e+3);
                break;
            end
        end 
    end
    disp('      all jobs completed ');
    disp('      time elapse= '+string(etime(clock(),c2))+' seconds');

    // collect all data ------------------------------------------------
    disp('  => collecting calculated results')
    kif_sbn.k_band=zeros(length(kif_sbn.SelBan),length(kif.k_point(:,1)));
    r1=1;
    for n=1:kif_sbn.Thread
        load(work_dir+'kif_sbn/k_band_'+string(n)+'.sod');
        r2=r1+length(k_band(1,:))-1;
        kif_sbn.k_band(:,r1:r2)=k_band;
        r1=r2+1;
        clear k_band;
    end

    // delete temporary folder -----------------------------------------
    disp('  => removing scratch folder')
    removedir(work_dir+'kif_sbn');

    // output information ==============================================
    disp('{kif_sbn} output information ...');
    fid(1)=mopen(project_name+'_kif_sbn.plb','a+');
    PIL_print_mat('kif_sbn.k_mesh, @f:f, 2D BZ k-mesh', kif_sbn.k_mesh,..
    'i',fid(1),'off'); 
    if kif_sbn.Verbosity=='more' then
        PIL_print_mat('kif_sbn.k_band, @f:f, selected eigenvalues,'+..
        ' [E(k1),E(k2),...]', kif_sbn.k_band,'r',fid(1));
    end
    mclose(fid(1));

    // finishing program ===============================================
    save(project_name+'_kif_sbn.sod','kif_sbn');
    disp('{kif_sbn}: finishing {kif_sbn} calculation ...');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction
