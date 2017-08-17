// **** Purpose ****
// PiLab density of state solver (level 3)
// **** variables ****
// << PiLab inputs >>
// [dsa.Mesh]: 1x3/1x2/1x1, int
// <= k-space mesh, the denser the better, but cost time
// [dsa.Interval]: 1x1, real
// <= energy interval, too small may cause delta functions
// [dsa.Draw]: 1x1, string, 'on' / 'off'
// <= whether draw DOS
// [dsa.Shift]: 1x1, string, 'on' / 'off'
// whether shift Ef=0 in plot
// << PiLab outputs >>
// [dsa.E_dsa]: tot energy division x 2, real
// => [E_point,normalized E_DOS]  
// **** Version ****
// 05/26/2014 1st version
// 05/12/2015 change reload process
// **** Comment ****


function PiLab_dsa(project_name)
    tic();
    disp('Starting {dsa} calculation ...');
    disp('=========== Message ===========');

    // loading variables ===============================================
    disp('## loading variables ...');
    PiLab_loader(project_name,'dsa','user','trim');
    load(project_name+'_dsa.sod'); 
    load(project_name+'_lat.sod');
    load(project_name+'_hop.sod');
    load(project_name+'_scc.sod');

    // variable check ===================================================
    disp('## checking variables ...')
    check_var=(length(dsa.Mesh)==length(lat.Sublatt(1,:))...
    & sum(dsa.Mesh==0)==0)
    if check_var~=%t then
        disp('Error: PiLab_dsa, dsa.Mesh has wrong dimension or 0 mesh!');
        abort;
    end
    check_var=(dsa.Draw=='on' | dsa.Draw=='off')
    if check_var~=%t then
        disp('Error: PiLab_dsa, dsa.Draw can be ''on'' or ''off'' !');
        abort;
    end
    check_var=(dsa.Shift=='on' | dsa.Shift=='off')
    if check_var~=%t then
        disp('Error: PiLab_dsa, dsa.Draw can be ''on'' or ''off'' !');
        abort;
    end
    // Core Part =======================================================
    disp('## running core part ...');
    // calculate DOS
    k_point=PIL_k_mesh(lat.recip_vec,dsa.Mesh);
    tot_state=length(hop.state_info(:,1));
    tot_k=length(k_point(:,1));
    E_level=zeros(tot_state,tot_k);
    for n=1:tot_k
        Hk=PIL_Hk_gen(k_point(n,:),lat.surr_site,...
        hop.state_info,scc.H_onsite,hop.hop_mat);
        E_level(:,n)=spec(Hk);
    end
    [E_point,E_DOS]=PIL_DOS(E_level(:),dsa.Interval);
    dsa.E_dsa=cat(2,E_point',E_DOS');

    // plot dos
    if dsa.Draw=='on' then
        xdel(winsid());
        if dsa.Shift=='on' then
            Ef=scc.E_Fermi;
        else
            Ef=0;    
        end
        plot(dsa.E_dsa(:,1)-Ef,dsa.E_dsa(:,2));
        plot([linspace(Ef,Ef,10)]',[linspace(0,1,10)]','r:');
        title('DOS',"fontsize", 3);
        xlabel('Energy',"fontsize", 3); 
        ylabel('DOS (arbitary unit)',"fontsize", 3);

    end

    // output information ==============================================
    disp('## output information ...')
    if dsa.Draw=='on' then
        xsave(project_name+'_dsa.scg');
        disp('   Output plot has been saved to '+project_name+'_dsa.scg'); 
    end

    // write to file
    fid(1)=mopen(project_name+'_dsa.plb','a+');
    PIL_print_mat('dsa.E_dsa, @full, [E_point,E_DOS]',dsa.E_dsa,'r',fid(1));
    mclose(fid(1));

    // finishing program ===============================================
    save(project_name+'_dsa.sod','dsa');
    disp('=============================');
    disp('Finishing {dsa} calculation ...');
    disp('# time elapse= '+string(toc())+ ' seconds');
endfunction
