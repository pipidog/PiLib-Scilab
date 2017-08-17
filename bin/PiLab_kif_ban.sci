// **** Purpose ****
// {kif_ban}: band structure along a particular k-path
// **** Variables ****
// ==== << PiLab inputs >> ====
// [kif_ban.SelBan]: 1xn, integer
// <= band index of the selected bands
// [kif_ban.StateProj]: 1xn, integer
// <= projection onto assigined states. use negative as divider
// [kif_ban.Verbosity]: 1x1, string, 'more' / 'less'
// <= whether to print out all eigenvalues and weighting 
// ==== << PiLab outputs >> ====
// [kif_ban.k_div]: 1x3, integer
// => number of k-pt of each k-path
// [kif_ban.k_band]: tot_SelBan * tot_k
// => band structure of each band structure
// [kif_ban.k_weight]: tot_SelBan * tot_k * tot_proj
// => projected weighting of each band states 
// **** Version ****
// 05/16/2016: first built
// **** Comment ****

function PiLab_kif_ban(project_name)
    disp('{kif_ban}: starting calculation ...');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{kif_ban}: loading variables ...');
    PiLab_loader(project_name,'kif_ban','user','trim');
    load(project_name+'_kif_ban.sod');
    load(project_name+'_kif.sod');
    disp('  => all variables loaded');

    // variable check ===================================================
    disp('{kif_ban}: checking variables ...');
    // check variables -----------------------------------------------------
    if kif.KptGen~='path' then
        disp('Error: PiLab_kif_ban, not working for kif.KptGen~=''path'' case');
        abort
    end
    if kif_ban.Verbosity~='more' & kif_ban.Verbosity~='less' then
        disp('Error: PiLab_kif_ban, kif_ban.Verbosity must be ''more'' or ''less''');
        abort
    end
    disp('  => all variables passed !');

    // core Part ========================================================
    disp('{kif_ban} running core part ...');
    work_dir=PIL_dir_path(pwd());

    // define necessary parameters
    load(work_dir+'kif/Ek/Ek_1.sod');
    tot_state=length(Ek(:,1));
    tot_k=length(kif.k_point(:,1));
    tot_sel=length(kif_ban.SelBan);
    tot_proj=length(find(kif_ban.StateProj<0));

    // built tasks of projected states
    if kif_ban.StateProj~=[] then
        if max(kif_ban.StateProj)> tot_state then
            disp('Error: PiLab_kif_ban, kif_ban.StateProj has non-exist stetes');
            abort
        end
        proj_ind=find(kif_ban.StateProj<0);
        state_proj=list();
        for n=1:length(proj_ind)
            r1=proj_ind(n)+1;
            if n~=length(proj_ind) then
                r2=proj_ind(n+1)-1;
            else
                r2=length(kif_ban.StateProj);
            end        
            state_proj(n)=kif_ban.StateProj(r1:r2);
        end
    end



    kif_ban.k_div=kif.k_div;
    kif_ban.k_band=zeros(tot_sel,tot_k);
    kif_ban.k_weight=zeros(tot_sel,tot_k,tot_proj);
    // calculate projected dos 
    disp('  => analyze each allocation\n');
    for n=1:kif.Allocation
        printf('     running %4d-th allocation\n',n);
        load(work_dir+'kif/Ek/Ek_'+string(n)+'.sod');
        kif_ban.k_band(:,kif.k_alloc(n,2):kif.k_alloc(n,3))=..
        Ek(kif_ban.SelBan,:);

        load(work_dir+'kif/Vk/Vk_'+string(n)+'.sod');
        if kif_ban.StateProj~=[] then
            // generate selected column index of Vk 
            idx=(repmat(kif_ban.SelBan,kif.k_alloc(n,1),1)..
            +repmat([0:kif.k_alloc(n,1)-1]'*tot_state,1,tot_sel))';
            idx=idx(:);

            for m=1:tot_proj
                Vk=abs(Vk(state_proj(m),idx)).^2
                den=zeros(1,length(Vk(:,2)));
                for p=1:length(Vk(:,1))
                    den=den+Vk(p,:)
                end
                kif_ban.k_weight(:,kif.k_alloc(n,2):kif.k_alloc(n,3),m)..
                =matrix(den,tot_sel,kif.k_alloc(n,1));
            end
        end
    end
    if kif_ban.StateProj==[] then
        proj_ind=1;
        ban.k_weight=ones(tot_sel,tot_k,1);
    end
    clear Ek Vk;

    // output information ==============================================
    disp('{kif_ban} output information ...');
    fid(1)=mopen(project_name+'_kif_ban.plb','a+');
    PIL_print_mat('kif_ban.k_div, @f:f, k_div', kif_ban.k_div..
    ,'i',fid(1),'off');
    if kif_ban.Verbosity=='more' then
        PIL_print_mat('kif_ban.k_band, @f:f, slected bands [En(k1), En(k2), etc.]'..
        ,kif_ban.k_band,'r',fid(1),'off');
        for n=1:tot_proj
            PIL_print_mat('kif_ban.k_weght(:,:,n), @f:f, '+string(n)..
            +'-th band-state-weighting',squeeze(kif_ban.k_weight(:,:,n))..
            ,'r',fid(1),'off');
        end
    end 
    mclose(fid(1));

    // finishing program ===============================================
    save(project_name+'_kif_ban.sod','kif_ban');
    disp('{kif_ban}: finishing {kif_ban} calculation ...');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction


