// This function collects all the core parts of PiLab. 
// It is convenient for users to implement PiLab functions in traditional
// Scilab executation runs. 
// Note: One should not change this file in any ways. All the content
//       of this file should be merely a copy and past of the core part
//       of each PiLab functions. 
function lat=PiLab_core_lat(lat)
    // generate reciprocal lattice vector
    [lat.recip_vec]=PIL_recip_vec(lat.Const*lat.Primitive);
    // generate lat.surr_site
    [lat.surr_site]=PIL_uc_nb(lat.Const*lat.Primitive...
    ,lat.Const*lat.Sublatt,lat.VecOrder,lat.NNCriterion,lat.NNOrder);
endfunction

function hop=PiLab_core_hop(lat,hop)
    // enforce orbitab identifier used in PIL_hop_mat
    hop.SiteOrb=hop.SiteOrb(:,[1,1,2]);  
    // to be the input order
    hop.SiteOrb(:,2)=[1:length(hop.SiteOrb(:,2))]';     

    hop.SiteOrb=gsort(hop.SiteOrb,'lr','i');
    [hop.state_info,hop.hop_mat]=PIL_hop_mat(lat.surr_site,hop.SiteOrb...
    ,hop.SKint,hop.Basis,hop.Order,hop.Filter);

    //generate LS coupling ---------------------------------
    hop.LS_mat=[];
    for n=1:length(hop.SiteOrb(:,1))
        select hop.SiteOrb(n,3)
        case 0
            hop.LS_mat=PIL_dirsum(hop.LS_mat...
            ,PIL_LS_coup('s',hop.LS,hop.Basis));
        case 1
            hop.LS_mat=PIL_dirsum(hop.LS_mat...
            ,PIL_LS_coup('p',hop.LS,hop.Basis));
        case 2
            hop.LS_mat=PIL_dirsum(hop.LS_mat...
            ,PIL_LS_coup('d',hop.LS,hop.Basis));
        case 3
            hop.LS_mat=PIL_dirsum(hop.LS_mat...
            ,PIL_LS_coup('f',hop.LS,hop.Basis));
        end
    end

    // generate  onsite energy ---------------
    hop.onsite_E=zeros(hop.LS_mat);

    // pick up selected sub-orbitals ---------
    if length(hop.SelState)~=0 then
        for n=1:size(hop.hop_mat)
            hop.hop_mat(n)=hop.hop_mat(n)(hop.SelState,hop.SelState,:);
        end           
        hop.LS_mat=hop.LS_mat(hop.SelState,hop.SelState);
        hop.onsite_E=diag(hop.OnsiteE);
        hop.state_info=hop.state_info(hop.SelState,:);
        // relabel new picked states
        hop.state_info(:,1)=(1:length(hop.SelState))';
    end
    // generate size of hop_mat
    hop.hop_size=zeros(size(hop.hop_mat),4);
    for n=1:size(hop.hop_mat)
        hop.hop_size(n,:)=[n,size(hop.hop_mat(n))];
    end
    // generate hop.state_info_text
    suborb_list=PIL_suborb_list(hop.Basis);
    hop.state_info_text=cat(2,string(hop.state_info(:,1:4)),...
    suborb_list(hop.state_info(:,$)));
endfunction

function scc=PiLab_core_scc(lat,hop,scc)
    // define necessary parameters
    tot_uc=prod(scc.Mesh);
    tot_uc_state=length(scc.Charge);
    tot_e=sum(scc.Charge)*tot_uc;

    // generate initial DM_uc_in
    if isfield(scc,'DM_out') & isfield(scc,'U_mat') then
        DM_in=scc.DM_out;
    else
        DM_in=diag(scc.Charge);
    end

    // generate all k-points 
    k_point=PIL_k_mesh(lat.recip_vec,scc.Mesh);
    k_point=PIL_vec_3d(k_point);

    // starting self-consistent calculation
    // HDD size: 1.5M ~ 1e+5 matrix elements ! 
    // buffer_size: the # of k-points to store in each file
    buffer_size=fix((scc.Memory/1.5)*1e+5/(tot_uc_state^2));
    if buffer_size > tot_uc
        buffer_size=tot_uc;
    end

    convergence=100;
    count=0;
    disp('   total eigenstate buffer files='+string(ceil(tot_uc/buffer_size))); 
    while convergence > scc.Converge
        count=count+1;
        disp('   < iteration='+string(count)+' >'); 
        // calculate U_mat
        scc.U_mat=zeros(tot_uc_state,tot_uc_state);
        for n=1:length(scc.HubU(:,1))
            scc.U_mat(scc.HubU(n,1),scc.HubU(n,1))...
            =scc.HubU(n,2)*(1/2-DM_in(scc.HubU(n,1),scc.HubU(n,1)));
        end
        // calculate Ef
        H_onsite=hop.onsite_E+scc.U_mat+hop.LS_mat;
        E_level=zeros(tot_uc_state,tot_uc);
        V_k=zeros(tot_uc_state,tot_uc_state,buffer_size);
        for n=1:tot_uc
            // diag Hk
            Hk=PIL_Hk_gen(k_point(n,:),lat.surr_site,hop.state_info...
            ,H_onsite,hop.hop_mat);
            [V,D]=spec(Hk);
            E_level(:,n)=real(diag(D));
            // save V to file
            n_buffer_mod=pmodulo(n,buffer_size);
            if n_buffer_mod==0 then
                V_k(:,:,$)=V;
                disp('   writing buffer file '+string(round(n/buffer_size)));
                save('V_k_'+string(round(n/buffer_size))+'.sod','V_k');
                V_k=zeros(V_k);
            elseif n==tot_uc 
                V_k(:,:,$)=V;
                disp('   writing buffer file '+string(ceil(n/buffer_size)));
                save('V_k_'+string(ceil(n/buffer_size))+'.sod','V_k');
            else
                V_k(:,:,n_buffer_mod)=V;
            end
        end
        clear V_k

        // calculate E_fermi
        scc.E_Fermi=PIL_Ef(E_level(:),tot_e,scc.Temperature);

        // calculate degenerate weighting
        state_weight=tot_e/length(find(E_level(:) <= scc.E_Fermi));

        //integrate all k-space to get new onsite density 
        scc.DM_out=zeros(tot_uc_state,tot_uc_state);
        for n=1:tot_uc
            tot_e_k=length(find(E_level(:,n) <= scc.E_Fermi))*state_weight;
            // read eigenvectors
            n_buffer_mod=pmodulo(n,buffer_size);
            select n_buffer_mod
            case 0
                V_tmp=V_k(:,:,$);
            case 1
                load('V_k_'+string(ceil(n/buffer_size))+'.sod');
                V_tmp=V_k(:,:,n_buffer_mod);
            else
                V_tmp=V_k(:,:,n_buffer_mod);
            end
            scc.DM_out=scc.DM_out...
            +PIL_DM_gen(E_level(:,n),0.01,scc.E_Fermi,tot_e_k,V_tmp)/tot_uc;
        end
        scc.DM_out=real(diag(diag(scc.DM_out)));
        tot_e_diff=abs(sum(diag(scc.DM_out))-sum(scc.Charge))/sum(scc.Charge)
        if tot_e_diff >= 10^-3 then
            disp('Warning: PiLab_scc, total charge difference='...
            +string(tot_e_diff*100)+'% > 0.1%');
        end


        // check convergence
        convergence=max(abs(scc.DM_out-DM_in));
        disp('   => Ef='+string(scc.E_Fermi)+' ; state_weight='...
        +string(state_weight)+' ; convergence='+string(convergence));
        DM_in=(1-scc.Mixing)*DM_in+scc.Mixing*scc.DM_out;
    end
    // detete temp files
    disp('   All eigenstate temporary files erased')
    for n=1:ceil(tot_uc/buffer_size)
        mdelete('V_k_'+string(n)+'.sod');
    end

    // write SCF information
    if convergence > scc.Converge then
        disp('   Exceed Maximal Iterations!');
    else
        disp('   Self-consistency reached!');
    end

    // Band Gap & DOS calculation
    E_level=gsort(E_level(:),'g','i');
    Ef_level=max(find(E_level <= scc.E_Fermi));
    scc.E_gap=E_level(Ef_level+1)-E_level(Ef_level);
    disp('Band Gap='+string(scc.E_gap));

    // define H_onsite
    scc.H_onsite=hop.onsite_E+scc.U_mat+hop.LS_mat;    
endfunction

function ban=PiLab_core_ban(lat,hop,scc,ban)
    // generate k-points
    k_path=zeros(ban.Path);
    if ban.Format=='coefficient' then
        for n=1:length(ban.Path(:,1)) // total points
            for m=1:length(ban.Path(1,:)) // dimension
                k_path(n,:)=k_path(n,:)+ban.Path(n,m)*lat.recip_vec(m,:);
            end
        end
    else
        k_path=ban.Path;
    end
    k_path=PIL_vec_3d(k_path);

    // calculate band structure
    [ban.k_point,ban.k_path_div]...
    =PIL_k_path(k_path,ban.Div,ban.DivType,lat.Const);
    tot_state=length(hop.state_info(:,1));
    tot_k=length(ban.k_point(:,1));
    ban.k_band=zeros(tot_state,tot_k);
    ban.k_vec=zeros(tot_state,tot_state,tot_k);
    for n=1:tot_k;
        Hk=PIL_Hk_gen(ban.k_point(n,:),lat.surr_site,...
        hop.state_info,scc.H_onsite,hop.hop_mat,ban.Conv);
        [V,D]=spec(Hk);
        ban.k_band(:,n)=clean(real(diag(D)));
        ban.k_vec(:,:,n)=V;
    end 

    // draw band structure
    if ban.Draw=='on' then
        xdel(winsid());
        if ban.Shift=='on' then
            ban.k_band=ban.k_band-scc.E_Fermi;
            scc.E_Fermi=0;
        end
        plot([1:tot_k]',[ban.k_band]','b');
        plot([1:tot_k]',[linspace(scc.E_Fermi,scc.E_Fermi,tot_k)]','r:');
        if length(ban.Path(:,1)) > 2 then
            for n=1:length(ban.Path(:,1))-2
                ban_max=max(ban.k_band);
                ban_min=min(ban.k_band);
                ban_width=max(ban.k_band)-min(ban.k_band)
                plot(sum(ban.k_path_div(1:n))*ones(10,1),...
                [linspace(ban_min-0.05*ban_width,ban_max+0.05*ban_width,10)]','k:');
            end
        end
        if length(ban.Ebound)==2
            a=gca();
            a.data_bounds=[1, ban.Ebound(1)...
            ;sum(ban.k_path_div), ban.Ebound(2)];
        end
        title('Band Structure','fontsize',4); 
        xlabel('$k$','fontsize',4); ylabel('Energy',"fontsize", 4); 
    end    
endfunction

function ban_den=PiLab_core_ban_den(hop,ban,ban_den)
    state_info=hop.state_info; 
    state_info_text=hop.state_info_text; 
    clear hop

    tot_k=length(ban.k_band(1,:));
    tot_check=length(ban_den.State(:,1));

    tot_suborb=length(state_info(:,1));
    tot_orb=max(state_info(:,3));
    tot_site=max(state_info(:,2));

    ban_den.site_den=zeros(tot_site,tot_check);
    ban_den.orb_den=zeros(tot_orb,tot_check);
    ban_den.suborb_den=zeros(tot_suborb,tot_check);
    ban_den.site_up_den=zeros(tot_site,tot_check);
    ban_den.site_dn_den=zeros(tot_site,tot_check);
    ban_den.orb_up_den=zeros(tot_orb,tot_check);
    ban_den.orb_dn_den=zeros(tot_orb,tot_check);


    for n=1:tot_check
        disp('   running '+string(n)+'-th band states')
        ban_den.suborb_den(:,n)=(abs(ban.k_vec(:,ban_den.State(n,2)...
        ,ban_den.State(n,1)))).^2;
        for m=1:tot_suborb
            site=state_info(m,2);
            orb=state_info(m,3);
            if grep(state_info_text(m,5),',u')==1 then
                ban_den.site_up_den(site,n)=ban_den.site_up_den(site,n)...
                +ban_den.suborb_den(m,n);
                ban_den.orb_up_den(orb,n)=ban_den.orb_up_den(orb,n)...
                +ban_den.suborb_den(m,n);            
            elseif grep(state_info_text(m,5),',d')==1
                ban_den.site_dn_den(site,n)=ban_den.site_dn_den(site,n)...
                +ban_den.suborb_den(m,n);
                ban_den.orb_dn_den(orb,n)=ban_den.orb_dn_den(orb,n)...
                +ban_den.suborb_den(m,n); 
            end
            ban_den.site_den(site,n)=ban_den.site_den(site,n)...
            +ban_den.suborb_den(m,n);
            ban_den.orb_den(orb,n)=ban_den.orb_den(orb,n)...
            +ban_den.suborb_den(m,n);
        end
        //check results
        if abs(sum(ban_den.suborb_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.suborb_den='...
            +string(sum(ban_den.suborb_den(:,n))));
        end
        if abs(sum(ban_den.orb_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.orb_den='...
            +string(sum(ban_den.orb_den(:,n))));
        end
        if abs(sum(ban_den.site_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.site_den='...
            +string(sum(ban_den.site_den(:,n))));
        end
        if abs(sum(ban_den.site_up_den(:,n)+ban_den.site_dn_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.site_spin='...
            +string(sum(ban_den.site_up_den(:,n)+ban_den.site_dn_den(:,n))))
        end
        if abs(sum(ban_den.orb_up_den(:,n)+ban_den.orb_dn_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.orb_spin='...
            +string(sum(ban_den.orb_up_den(:,n)+ban_den.orb_dn_den(:,n))));
        end
    end
    // Draw results  
    if ban_den.Draw=='on'
        for n=1:tot_check
            subplot(tot_check,1,n);
            select ban_den.DrawType
            case 'site'
                xlabel('site','fontsize',3);
                ylabel(string(n)+' den','fontsize',3);
                bar(ban_den.site_den(:,n));
            case 'orb'
                xlabel('orb','fontsize',3);
                ylabel(string(n)+' den','fontsize',3);
                bar(ban_den.orb_den(:,n));
            case 'suborb'
                xlabel('suborb','fontsize',3);
                ylabel(string(n)+' den','fontsize',3);
                bar(ban_den.suborb_den(:,n)); 
            case 'site_up'
                xlabel('site','fontsize',3);
                ylabel(string(n)+' up den','fontsize',3);
                bar(ban_den.site_up_den(:,n));
            case 'site_dn'
                xlabel('site','fontsize',3);
                ylabel(string(n)+' dn den','fontsize',3);
                bar(ban_den.site_dn_den(:,n));
            case 'orb_up'
                xlabel('orb','fontsize',3);
                ylabel(string(n)+' up den','fontsize',3);
                bar(ban_den.orb_up_den(:,n));
            case 'orb_dn'
                xlabel('orb','fontsize',3);
                ylabel(string(n)+' dn den','fontsize',3);
                bar(ban_den.orb_dn_den(:,n));
            end
            set(gcf(),'background',8)
        end
    end 
endfunction

function dsa=PiLab_core_dsa(lat,hop,scc,dsa)
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
endfunction

function zti=PiLab_core_zti(lat,hop,scc,zti)
    // check TR existence ----------------------------------------------
    disp('   Checking time-reversal symmetry');
    if max(abs(PIL_TR_check(lat,hop,scc))) >=1e-4
        disp('Error: PiLab_zti, time-reversal check does not pass!');
        abort;
    else
        disp('   Passed!');
    end

    // calculate Z2 ----------------------------------------------------    
    lat_size=size(lat.Primitive);
    tot_k_plane=prod(2*zti.Mesh);
    select lat_size(1)
    case 2
        [zti.Z2_val,zti.n_field]...
        =PIL_Z2_cal(lat,hop,scc,[],[3,0],zti.Mesh,zti.OccBand);
        disp('   Z2_val='+string(zti.Z2_val));

    case 3
        // Z2 calculation ----------------------------------------------
        zti.Z2_plane=zeros(3,2);
        zti.n_field=zeros(tot_k_plane,6,3,2);
        for n=1:3
            for m=0:1
                [Z2_plane,n_plane]...
                =PIL_Z2_cal(lat,hop,scc,[],[n,m],zti.Mesh,zti.OccBand);
                zti.Z2_plane(n,m+1)=Z2_plane;
                zti.n_field(:,:,n,m+1)=n_plane;
                disp('   Z2_val of plane ['+string(n)+','+string(m)+']='...
                +string(zti.Z2_plane(n,m+1)));
            end
        end
        // check Z2_val ------------------------------------------------
        Z2_check=pmodulo(zti.Z2_plane(:,1)+zti.Z2_plane(:,2),2);
        if sum(Z2_check)~=0 & sum(Z2_check)~=3
            disp('Error: PiLab_zti, sum of zti.Z2_plane inconsistent!');
            abort;
        end
        zti.Z2_val=pmodulo(...
        round([sum(zti.Z2_plane(1,:)),zti.Z2_plane(:,2)']),2);
        disp('   -----------------');
        disp('   total Z2_val= '+string(zti.Z2_val));
    end

    // plot n_field ------------------------------------------------
    select lat_size(1)
    case 2
        if zti.Draw=='on' then
            n_plot=round(matrix(zti.n_field(:,3),2*zti.Mesh(1)...
            ,2*zti.Mesh(2)));
            n_plot(n_plot==-1)=2; 
            n_plot(n_plot==0)=18; 
            n_plot(n_plot==1)=5;
            Matplot(n_plot);
            title('n_field','fontsize',4);

        end
    case 3
        if zti.Draw=='on'
            for n=1:3
                for m=1:2
                    n_plot=round(matrix(zti.n_field(:,3,n,m),...
                    2*zti.Mesh(1),2*zti.Mesh(2)));
                    n_plot(n_plot==-1)=2; 
                    n_plot(n_plot==0)=18; 
                    n_plot(n_plot==1)=5;
                    select m
                    case 1
                        subplot(2,3,n);
                    case 2
                        subplot(2,3,n+3);                    
                    end
                    Matplot(n_plot);
                    title('n_field @ ['+string(n)+','+string(m-1)+']'...
                    ,'fontsize',4);
                end
            end
        end
    end
endfunction

function chn=PiLab_core_chn(lat,hop,scc,chn)
    [chn.k_point,chn.Chern_val,chn.Fk_field]...
    =PIL_Chern_cal(lat,hop,scc,[],chn.Mesh,chn.OccBand,chn.Kdiff)
endfunction

function flq=PiLab_core_flq(lat,hop,scc,flq)
    // calculate renormalized matrix elements -------------------------
    // flq.H_onsite(order_index)
    // flq.hop_mat(order_index)(sublatt_index)
    flq.H_onsite=list()
    flq.hop_mat=list();
    for p=1:flq.Order+1 // run photon process order, 0~flq.Order
        // Oniste renormalization
        flq.H_onsite(p)=scc.H_onsite;
        flq.H_onsite(p)=flq.H_onsite(p)...
        *PIL_flqint(flq.Frequency,p-1,flq.Amplitude,flq.Phase,[0,0,0]);

        // hopping renormalization
        flq.hop_mat(p)=list();
        for n=1:length(hop.hop_size(:,1)) // run all sublattice
            r1=lat.surr_site(n)(1,7:9);
            flq.hop_mat(p)(n)=hop.hop_mat(n);
            for m=1:hop.hop_size(n,4) // run all nearest neighbor
                // renormalize hopping integrals
                r2=lat.surr_site(n)(m+1,7:9);    
                flq.hop_mat(p)(n)(:,:,m)=flq.hop_mat(p)(n)(:,:,m)...
                *PIL_flqint(flq.Frequency,p-1...
                ,flq.Amplitude,flq.Phase,r2-r1);
            end
        end
    end

    // generate Floquet state index ------------------------------------
    flq.state_info=[];
    for n=-flq.Order:flq.Order
        flq.state_info=cat(1,flq.state_info...
        ,cat(2,n*ones(hop.state_info(:,1)),hop.state_info(:,2:$)));
    end
    flq.state_info=cat(2,[1:length(flq.state_info(:,1))]'...
    ,flq.state_info);
    // state info in text
    // generate flq.state_info_text
    suborb_list=PIL_suborb_list(hop.Basis);
    flq.state_info_text=cat(2,string(flq.state_info(:,1:5)),...
    suborb_list(flq.state_info(:,$)))

    // generate size of flq.hop_mat for PiLab_loader -------------------
    flq.hop_size=[];
    for p=1:flq.Order+1
        flq.hop_size=cat(1,flq.hop_size,cat(2,...
        p*ones(length(hop.hop_size(:,1)),1),hop.hop_size));
    end

    // genreate Ef and Floquet band gap
    k_point=PIL_k_mesh(lat.recip_vec,flq.Mesh);
    k_point=PIL_vec_3d(k_point);
    tot_k=length(k_point(:,1));
    tot_state=length(flq.state_info(:,1));
    tot_e=flq.Charge*tot_k;

    k_band=zeros(tot_state,tot_k);
    for n=1:tot_k
        Hk=PIL_Hk_flq(k_point(n,:),lat,flq,'full');
        [V,D]=spec(Hk);
        k_band(:,n)=real(clean(diag(D)));
    end
    k_band=k_band(:);
    k_band=gsort(k_band,'g','i');
    flq.E_gap=k_band(tot_e+1)-k_band(tot_e)
    flq.E_Fermi=PIL_Ef(k_band,tot_e,flq.Temp);
endfunction

function flq_ban=PiLab_core_flq_ban(lat,flq,flq_ban)
    // generate k-paths if coefficient ---------------------------------
    k_path=zeros(flq_ban.Path);
    if flq_ban.Format=='coefficient' then
        for n=1:length(flq_ban.Path(:,1)) // total points
            for m=1:length(flq_ban.Path(1,:)) // dimension
                k_path(n,:)=k_path(n,:)...
                +flq_ban.Path(n,m)*lat.recip_vec(m,:);
            end
        end
    else
        k_path=flq_ban.Path;
    end
    k_path=PIL_vec_3d(k_path);

    // generate k-points -----------------------------------------------
    [flq_ban.k_point,flq_ban.k_path_div]...
    =PIL_k_path(k_path,flq_ban.Div,flq_ban.DivType,lat.Const);
    tot_state=length(flq.state_info(:,1));
    tot_k=length(flq_ban.k_point(:,1));
    flq_ban.k_band=zeros(tot_state,tot_k);
    //flq_ban.k_vec=zeros(tot_state,tot_state,tot_k);

    // calculate band structure ----------------------------------------
    for n=1:tot_k
        if fix(10*n/tot_k)~=fix(10*(n-1)/tot_k) then
            disp('   Floquet Band Calculation Process= '...
            +string(round(100*n/tot_k))+'%');
        end
        Hk=PIL_Hk_flq(flq_ban.k_point(n,:),lat,flq,flq_ban.Conv);
        [V,D]=spec(Hk);
        flq_ban.k_band(:,n)=real(clean(diag(D)));
        //flq_ban.k_vec(:,:,n)=clean(V)
        // print k_vec, so don't need to store them in memory
    end
    // draw band structure
    if flq_ban.Draw=='on' then
        xdel(winsid());
        if flq_ban.Shift=='on' then
            flq_ban.k_band=flq_ban.k_band-flq.E_Fermi;
            flq.E_Fermi=0;
        end
        plot([1:tot_k]',[flq_ban.k_band]','b');
        plot([1:tot_k]',[linspace(flq.E_Fermi,flq.E_Fermi,tot_k)]','r:');
        if length(flq_ban.Path(:,1)) > 2 then
            for n=1:length(flq_ban.Path(:,1))-2
                ban_max=max(flq_ban.k_band);
                ban_min=min(flq_ban.k_band);
                ban_width=max(flq_ban.k_band)-min(flq_ban.k_band)
                plot(sum(flq_ban.k_path_div(1:n))*ones(10,1),...
                [linspace(ban_min-0.05*ban_width...
                ,ban_max+0.05*ban_width,10)]','k:');
            end
        end
        if length(flq_ban.Ebound)==2
            a=gca();
            a.data_bounds=[1, flq_ban.Ebound(1)...
            ;sum(flq_ban.k_path_div), flq_ban.Ebound(2)];
        end
        title('Floquet Band Structure','fontsize',4); 
        xlabel('$k$','fontsize',4); ylabel('Energy',"fontsize", 4);
    end
endfunction

function flq_zti=PiLab_core_flq_zti(lat,hop,flq,flq_zti)
    disp('## running core part ...');
    // check TR existence 
    disp('   Checking time-reversal symmetry...');
    if max(abs(PIL_TR_check(lat,flq))) >= 1e-4
        disp('Error: PiLab_flq_zti, time-reversal check does not pass!');
        abort;
    else
        disp('   Passed!');
    end

    // run Z2 calculation
    lat_size=size(lat.Primitive);
    tot_k_plane=prod(2*flq_zti.Mesh);
    select lat_size(1)
    case 2
        [flq_zti.Z2_val,flq_zti.n_field]...
        =PIL_Z2_cal(lat,[],[],flq,[3,0],flq_zti.Mesh,flq_zti.OccBand);
        disp('   Z2_val='+string(flq_zti.Z2_val));
    case 3
        // Z2 calculation ----------------------------------------------
        flq_zti.Z2_plane=zeros(3,2);
        flq_zti.n_field=zeros(tot_k_plane,3,3,2);
        for n=1:3
            for m=0:1
                [Z2_plane,n_plane]...
                =PIL_Z2_cal(lat,[],[],flq,[n,m],flq_zti.Mesh...
                ,flq_zti.OccBand);
                flq_zti.Z2_plane(n,m+1)=Z2_plane;
                flq_zti.n_field(:,:,n,m+1)=n_plane;
                disp('   Z2_val of plane ['+string(n)+','+string(m)+']='...
                +string(flq_zti.Z2_plane(n,m+1)));
            end
        end
        // check Z2_val ------------------------------------------------
        Z2_check=pmodulo(flq_zti.Z2_plane(:,1)+flq_zti.Z2_plane(:,2),2);
        if sum(Z2_check)~=0 & sum(Z2_check)~=3
            disp('Error: PiLab_flq_zti, sum of flq_zti.Z2_plane inconsistent!');
            abort;
        end
        flq_zti.Z2_val=pmodulo(...
        round([sum(flq_zti.Z2_plane(1,:)),flq_zti.Z2_plane(:,2)']),2);
        disp('   -----------------');
        disp(['   total Z2_val=',strcat(string(flq_zti.Z2_val))]);
    end
endfunction

function flq_chn=PiLab_core_flq_chn(lat,flq,flq_chn)
    [flq_chn.k_point,flq_chn.Chern_val,flq_chn.Fk_field]...
    =PIL_Chern_cal(lat,[],[],flq,flq_chn.Mesh,flq_chn.OccBand,flq_chn.Kdiff);
endfunction

function flq_ban_den=PiLab_core_flq_ban_den(flq,flq_ban,flq_ban_den)
    state_info=flq.state_info; 
    state_info_text=flq.state_info_text; 
    if grep(flq.state_info_text(1,$),',u')==1 ...
        | grep(flq.state_info_text(1,$),',d')==1 then
        basis='spin'
    else
        basis=' '
    end
    flq_order=flq.Order;
    clear flq

    tot_k=length(flq_ban.k_band(1,:));
    tot_check=length(flq_ban_den.State(:,1));

    tot_order=2*flq_order+1;
    tot_flq_suborb=length(state_info(:,1));
    tot_suborb=tot_flq_suborb/tot_order;
    tot_orb=max(state_info(:,4));
    tot_site=max(state_info(:,3));

    flq_ban_den.site_den=zeros(tot_site,tot_check);
    flq_ban_den.site_up_den=zeros(tot_site,tot_check);
    flq_ban_den.site_dn_den=zeros(tot_site,tot_check);
    flq_ban_den.orb_den=zeros(tot_orb,tot_order,tot_check);
    flq_ban_den.orb_up_den=zeros(tot_orb,tot_order,tot_check);
    flq_ban_den.orb_dn_den=zeros(tot_orb,tot_order,tot_check);
    flq_ban_den.suborb_den=zeros(tot_suborb,tot_order,tot_check);
    flq_ban_den.flq_den=zeros(tot_flq_suborb,tot_check);
    for n=1:tot_check
        disp('   running '+string(n)+'-th band states')
        flq_ban_den.flq_den(:,n)=(abs(...
        flq_ban.k_vec(:,flq_ban_den.State(n,2),flq_ban_den.State(n,1)))).^2;
        for m=1:tot_flq_suborb
            order=state_info(m,2)+flq_order+1;
            site=state_info(m,3);
            orb=state_info(m,4);
            suborb=pmodulo(m,tot_suborb)...
            +PIL_kdelta(pmodulo(m,tot_suborb),0)*tot_suborb;

            if grep(state_info_text(m,6),',u')==1 then
                flq_ban_den.site_up_den(site,n)=...
                flq_ban_den.site_up_den(site,n)...
                +flq_ban_den.flq_den(m,n);

                flq_ban_den.orb_up_den(orb,order,n)=...
                flq_ban_den.orb_up_den(orb,order,n)...
                +flq_ban_den.flq_den(m,n);            
            elseif grep(state_info_text(m,6),',d')==1
                flq_ban_den.site_dn_den(site,n)=...
                flq_ban_den.site_dn_den(site,n)...
                +flq_ban_den.flq_den(m,n);

                flq_ban_den.orb_dn_den(orb,order,n)=...
                flq_ban_den.orb_dn_den(orb,order,n)...
                +flq_ban_den.flq_den(m,n); 
            end
            flq_ban_den.site_den(site,n)=...
            flq_ban_den.site_den(site,n)...
            +flq_ban_den.flq_den(m,n);

            flq_ban_den.orb_den(orb,order,n)=...
            flq_ban_den.orb_den(orb,order,n)...
            +flq_ban_den.flq_den(m,n);

            flq_ban_den.suborb_den(suborb,order,n)=...
            flq_ban_den.suborb_den(suborb,order,n)...
            +flq_ban_den.flq_den(m,n);
        end
    end
    // Draw results  
    if flq_ban_den.Draw=='on'
        for n=1:tot_check
            select flq_ban_den.DrawType
            case 'site'
                subplot(tot_check,1,n);
                xlabel('site','fontsize',3);
                ylabel(string(n)+' den','fontsize',3);
                bar(flq_ban_den.site_den(:,n));
                set(gcf(),'background',8)
            case 'site_up'
                subplot(tot_check,1,n);
                xlabel('site','fontsize',3);
                ylabel(string(n)+' up-den','fontsize',3);
                bar(flq_ban_den.site_up_den(:,n));
                set(gcf(),'background',8)
            case 'site_dn'
                subplot(tot_check,1,n);
                xlabel('site','fontsize',3);
                ylabel(string(n)+' dn-den','fontsize',3);
                bar(flq_ban_den.site_dn_den(:,n));
                set(gcf(),'background',8)
            case 'orb'
                for m=1:tot_order
                    subplot(tot_check,tot_order,(n-1)*tot_order+m);
                    bar(flq_ban_den.orb_den(:,m,n));
                    xlabel('orb ('+string(m-flq_order-1)+')','fontsize',3);
                    ylabel(string(n)+' den','fontsize',3);
                    set(gcf(),'background',8)
                end
            case 'orb_up'
                for m=1:tot_order
                    subplot(tot_check,tot_order,(n-1)*tot_order+m);
                    bar(flq_ban_den.orb_up_den(:,m,n));
                    xlabel('orb ('+string(m-flq_order-1)+')','fontsize',3);
                    ylabel(string(n)+' up-den','fontsize',3);
                    set(gcf(),'background',8)
                end
            case 'orb_dn'
                for m=1:tot_order
                    subplot(tot_check,tot_order,(n-1)*tot_order+m);
                    bar(flq_ban_den.orb_dn_den(:,m,n));
                    xlabel('orb ('+string(m-flq_order-1)+')','fontsize',3);
                    ylabel(string(n)+' dn-den','fontsize',3);
                    set(gcf(),'background',8)
                end
            case 'suborb'
                for m=1:tot_order
                    subplot(tot_check,tot_order,(n-1)*tot_order+m);
                    bar(flq_ban_den.suborb_den(:,m,n));
                    xlabel('suborb ('+string(m-flq_order-1)+')','fontsize',3);
                    ylabel(string(n)+' den','fontsize',3);
                    set(gcf(),'background',8)
                end
            end
        end
    end 
endfunction

