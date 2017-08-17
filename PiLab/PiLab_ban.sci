// **** Purpose ****
// PiLab band structure calculator (level 3)
//==== << PiLab inputs >> ====
//[ban.Format]: 1x1, string, 'red' / 'cart'
//<= format of ban.Path points, reduced or cartisian k-points
//[ban.Path]: nx3, n>=2, real            
//<= points to defined your paths, ban.Format defines their meaning 
//[ban.DivType]: 1x1, string, 'unit' / 'seg'
//<= how to divide k_path, 'unit': division for unit length
//'seg': divison for each segment  
//[ban.Div]: 1x1, int
//<= k points mesh of each path, ban.DivType defines their meaning
//[ban.Fermi]: 1x1, real
//<= value of Fermi energy, should be obtained from scc
//or ab initio if wan.
//[ban.Draw]: 1x1, string, 'on' / 'off'
//<= whether draw band structure
//[ban.Shift]: 1x1, string, 'on' / 'off'
//<= whether to shift Ef to 0 in plot
//[ban.Ebound]: 1x2, real, optional, default:[], i.e. natural range
//<= energy window to print band, [Emin, Emax], if ban.shift='on', 
//then Ebound assumes Ef=0.
//[ban.Verbosity]: 'less' / 'more'
//<= whether to print H(k) and eigenvectors
//
//==== << PiLab output >> ====
//[ban.k_path_div]: nx1, int
//=> number of divisions of each segment
//[ban.k_point]: tot_k x 4, real
//=> [kx,ky,kz]
//[ban.k_band]: tot_state x tot_k, real
//=> band energies, [En(k1),En(k2)...]
//[ban.H_k]: tot_k x tot_state x tot_state
//=> Hamiltonian matrix at each k-point
//[k_vec]: tot_k x tot_state x tot_state 
//=> eigenvectors of each k-point
// **** Version ****
// 05/29/2014 1st version
// 05/12/2015 change reload process
// 02/17/2016 change to fit PiLab 1.0 version
// **** Comment ****
// 1. This code generates:
//  1).the k-point scanned (ban.k_point(:,:))
//  2).the band energy (ban.k_band(:,:))
//  3).the eigenvector matrix at each k-point (k_vec(:,:,n))

function PiLab_ban(project_name)
    disp('{ban}: starting calculation ');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{ban}: loading variables ');
    PiLab_loader(project_name,'ban','user','trim');
    load(project_name+'_ban.sod');
    load(project_name+'_ham.sod');
    disp('  => all variables loaded')
    // variable check ===================================================
    disp('{ban}: checking variables ')
    check_var=(ban.Format=='red' | ban.Format=='cart');
    if check_var~=%t then
        disp('Error: PiLab_ban, ban.Format must be '...
        +'''red or ''cart'' !');
        abort;
    end
    check_var=(length(ban.Path(1,:))==3 & length(ban.Path(:,1)) >=2);
    if check_var~=%t then
        disp('Error: PiLab_ban, ban.Path has wrong input!');
        abort;
    end
    check_var=((ban.DivType=='unit') | (ban.DivType=='seg'));
    if check_var~=%t then
        disp('Error: PiLab_ban, ban.DivType must be ''unit'' or ''seg''');
        abort;
    end
    check_var=((ban.Draw=='on') | (ban.Draw=='off'));
    if check_var~=%t then
        disp('Error: PiLab_ban, ban.Draw has wrong input!');
        abort;
    end
    check_var=((ban.Shift=='on') | (ban.Shift=='off'));
    if check_var~=%t then
        disp('Error: PiLab_ban, ban.Shift has wrong input!');
        abort;
    end
    check_var= length(ban.Ebound)==0 | length(ban.Ebound)==2 
    if check_var~=%t then
        disp('Warning: PiLab_ban, ban.Ebound must [Emin, Emax] or [], '..
        +'set to default!');
        ban.Ebound=[]
    end
    disp('  => all variables passed')

    // core Part ========================================================
    disp('{ban} running core part ');
    // generate k-points
    disp('  => generating k-points')
    k_path=zeros(ban.Path);
    if ban.Format=='red' then
        [ban.k_point,ban.k_path_div]=..
        PIL_k_path(ban.Path*ham.rec_vec,ban.Div,ban.DivType);
    else
        [ban.k_point,ban.k_path_div]=..
        PIL_k_path(ban.Path,ban.Div,ban.DivType);
    end

    // prepare projection parameters 
    tot_state=length(ham.state_info(:,1));
    tot_k=length(ban.k_point(:,1));
    if ban.StateProj~=[] then
        if max(ban.StateProj)> tot_state then
            disp('Error: PiLab_ban, ban.StateProj has non-exist stetes');
            abort
        end
        // built tasks of projected states
        proj_ind=find(ban.StateProj<0);
        tot_proj=length(proj_ind);
        state_proj=list();
        for n=1:length(proj_ind)
            r1=proj_ind(n)+1;
            if n~=length(proj_ind) then
                r2=proj_ind(n+1)-1;
            else
                r2=length(ban.StateProj);
            end        
            state_proj(n)=ban.StateProj(r1:r2);
        end
        ban.k_weight=zeros(tot_state,tot_k,length(proj_ind));
    else
        proj_ind=1;
        ban.k_weight=ones(tot_state,tot_k,1);
        tot_proj=0;
    end

    // calculate band structure    
    disp('  => diagonalizing Hamiltonian')
    printf('\n')
    printf('     total k-points= %d\n',tot_k);
    printf('     total projection= %d\n',tot_proj);
    ban.k_band=zeros(tot_state,tot_k);
    c2=clock();
    for n=1:tot_k
        [hk,err]=PIL_Hk_R(ban.k_point(n,:),ham.lat_vec,ham.uc_index,[],ham.Hr_mat);
        [Vk,Dk]=spec(hk);
        ban.k_band(:,n)=diag(Dk);
        if ban.StateProj~=[] then
            Vk=abs(Vk).^2
            for m=1:tot_proj
                ban.k_weight(:,n,m)=sum(Vk(state_proj(m),:),'r')';
            end
        end
        clear Vk Dk hk;
        if tot_k>=10 & pmodulo(n,fix(tot_k/10))==0 then
            printf('     %4d k-points calculated, time=%f\n',..
            n,etime(clock(),c2));
            c2=clock();
        end
    end

    // shift band energy
    if ban.Shift=='on' then
        ban.k_band=ban.k_band-ban.Fermi;
        ban.Fermi=0;
    end 

    // draw band structure
    if ban.Draw=='on' then
        // draw normal band
        xdel(winsid());
        figure(1);
        a=gcf();
        a.background=-2;
        plot([1:tot_k]',[ban.k_band]','thickness',ban.Thickness,'color','b');
        plot([1:tot_k]',[linspace(ban.Fermi,ban.Fermi,tot_k)]','r:');
        if length(ban.Path(:,1)) > 2 then
            for n=1:length(ban.Path(:,1))-2
                ban_max=max(ban.k_band);
                ban_min=min(ban.k_band);
                ban_width=max(ban.k_band)-min(ban.k_band)
                plot(sum(ban.k_path_div(1:n))*ones(10,1),..
                [linspace(ban_min-0.05*ban_width,..
                ban_max+0.05*ban_width,10)]','k:');
            end
        end
        if length(ban.Ebound)==2
            a=gca();
            a.data_bounds=[1, ban.Ebound(1)...
            ;sum(ban.k_path_div), ban.Ebound(2)];
        end
        title('Band Structure','fontsize',4); 
        xlabel('$k$','fontsize',4); ylabel('Energy',"fontsize", 4);
        a=gca(); 
        a.tight_limits='on';
        a.thickness=3;
        a.font_size=3

        // draw projected bands
        if ban.StateProj~=[] then
            for n=1:length(proj_ind);
                figure(n+1);
                a=gcf();
                a.background=-2;

                // plot weighted band
                PIL_scatter_plot(matrix(repmat([1:tot_k],tot_state,1),-1,1),..
                ban.k_band(:),..
                matrix(ban.k_weight(:,:,n),-1,1),..
                ban.MarkerSize,flipdim(oceancolormap(64),1));
                                
                // plot divider and fermi level 
                // **scilab use stack, so childern(1) is always the first item
                plot([1:tot_k]',[linspace(ban.Fermi,ban.Fermi,tot_k)]','r:');  
                a=gca();
                a.children(1).children.thickness=1 

                if length(ban.Path(:,1)) > 2 then
                    ban_max=max(ban.k_band);
                    ban_min=min(ban.k_band);
                    ban_width=max(ban.k_band)-min(ban.k_band)
                    for m=1:length(ban.Path(:,1))-2
                        plot(sum(ban.k_path_div(1:m))*ones(10,1),..
                        [linspace(ban_min-0.05*ban_width,..
                        ban_max+0.05*ban_width,10)]','k:');
                        a.children(1).children.thickness=1

                    end
                end

                // set plot range
                if length(ban.Ebound)==2
                    a=gca();
                    a.data_bounds=[1, ban.Ebound(1)...
                    ;sum(ban.k_path_div), ban.Ebound(2)];
                end

                // set label
                title('Projected Band-'+string(n),'fontsize',4); 
                xlabel('$k$','fontsize',4); ylabel('Energy',"fontsize", 4);
            end
        end
    end

    // output information ==============================================
    disp('{ban} output information ')
    if ban.Draw=='on' then
        scf(1)
        xsave(project_name+'_ban.scg');
        disp('  => Output plot '+project_name+'_ban.scg saved'); 
        if ban.StateProj~=[] then
            for n=1:length(proj_ind)
                scf(n+1)
                xsave(project_name+'_ban_p'+string(n)+'.scg');
                disp('     Output plot '+project_name..
                +'_ban_p'+string(n)+'.scg saved');
            end
        end
    end
    fid(1)=mopen(project_name+'_ban.plb','a+');
    PIL_print_mat('ban.k_path_div, @f:f, division of each path'...
    , ban.k_path_div,'i',fid(1));
    PIL_print_mat('ban.k_point, @f:f, [kx,ky,kz]'...
    ,ban.k_point,'r',fid(1));
    PIL_print_mat('ban.k_band, @f:f, [En(k1),En(k2)...]'...
    ,ban.k_band,'r',fid(1));
    for n=1:length(proj_ind)
        PIL_print_mat('ban.k_weight(:,:,'+string(n)+'), @f:f, '..
        +string(n)+'-th projected weight of each eigenvalue'..
        ,ban.k_weight(:,:,n),'r',fid(1));
    end
    mclose(fid(1))
    // finishing program ===============================================
    save(project_name+'_ban.sod','ban');
    disp('{ban}: finishing {ban} calculation ');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction
