// **** Purpose ****
// load the variables of a PiLab file and save it to a .sod file
// it can also erase the non-user defined variables in PiLab files.
// **** Variables ****
// [project_name]: 1x1, string
// <= the project name
// [job_type]: 1x1, string
// <= job_type, ex: lat, hop, scc, ... 
// [read_task]: 1x1, string, 'user', 'all'
// <= 'user': read user defined variables
//    'all': read all variables
// [trim_task]: 1x1, string, 'keep', 'trim'
// <= 'keep': keep plb file unchange
//    'trim': erase calculated PiLib variables 
// **** Version ****
// 5/26/2014 first built
// 8/26/2014 add round function to lat.surr_site 
// **** Comment ****

function PiLab_loader(project_name,job_type,read_task,trim_task)
    [lhs,rhs]=argn();
    select rhs
    case 2
        read_task='all';
        trim_task='keep';
    case 3
        trim_task='keep';
    end

    if (read_task~='all') & (read_task~='user') then
        disp('Error: PiLab_loader, read_task input is wrong!');
        abort
    end
    if (trim_task~='keep') & (trim_task~='trim') then
        disp('Error: PiLab_loader, trim_task input is wrong!');
        abort
    end

    PiLab_data=0;
    fid=zeros(1,2);
    fid(1)=mopen(project_name+'_'+string(job_type)+'.plb','r');
    fid(2)=mopen(project_name+'_'+string(job_type)+'.tmp','w');

    // read user inputs ##############################################
    count_line=0;
    PiLib_variable='no';
    while meof(fid(1))==0
        read_data=mgetl(fid(1),1);
        if length(grep(read_data,'PiLib Variable'))==0 then
            count_line=count_line+1;
            if length(read_data)~=0 then
                mputl(read_data,fid(2));
            end
        else
            PiLib_variable='yes';
            break;
        end
    end
    mclose(fid(2));
    exec(project_name+'_'+job_type+'.tmp',-1);

    // read PiLib data ###############################################

    if (PiLib_variable=='yes') & (read_task=='all') then
        mseek(0,fid(1));
        mgetl(fid(1),count_line);
        select job_type
            // Hr construct ============================================
        case 'lat'
            lat.recip_vec=PIL_read_mat(fid(1));
            lat.surr_site=list();
            for n=1:length(lat.Sublatt(:,1))
                lat.surr_site(n)=PIL_read_mat(fid(1));
                lat.surr_site(n)(:,1)=round(lat.surr_site(n)(:,1));
                lat.surr_site(n)(:,3:6)=round(lat.surr_site(n)(:,3:6));
            end
        case 'hop'
            hop.state_info_text=PIL_read_mat(fid(1));
            hop.state_info=PIL_read_mat(fid(1));
            tmp=PIL_read_mat(fid(1));
            hop.LS_mat=PIL_sparse(tmp,'sparse','tri');
            tmp=PIL_read_mat(fid(1));
            hop.onsite_E=PIL_sparse(tmp,'sparse','tri');
            hop.hop_size=PIL_read_mat(fid(1));
            hop.hop_mat=list();
            for n=1:length(hop.hop_size(:,1))
                hop.hop_mat(n)=zeros(hop.hop_size(n,2),...
                hop.hop_size(n,3),hop.hop_size(n,4));
                for m=1:hop.hop_size(n,4)
                    tmp=PIL_read_mat(fid(1));
                    hop.hop_mat(n)(:,:,m)=PIL_sparse(tmp,'sparse','all');
                end
            end
        case 'scc'
            scc.E_Fermi=PIL_read_mat(fid(1));
            scc.E_gap=PIL_read_mat(fid(1));
            scc.DM_out=full(PIL_read_mat(fid(1)));
            scc.DM_out=scc.DM_out+scc.DM_out'-diag(diag(scc.DM_out));
            scc.U_mat=full(PIL_read_mat(fid(1)));
            scc.U_mat=scc.U_mat+scc.U_mat'-diag(diag(scc.U_mat));
        case 'flq'
            flq.state_info_text=PIL_read_mat(fid(1));
            flq.state_info=PIL_read_mat(fid(1));
            flq.E_Fermi=PIL_read_mat(fid(1));
            flq.E_gap=PIL_read_mat(fid(1));
            flq.H_onsite=list();
            flq.hop_mat=list();
            for p=1:flq.Order+1
                flq.H_onsite(p)...
                =PIL_sparse(PIL_read_mat(fid(1)),'sparse','tri');
                flq.hop_mat(p)=list();
            end
            flq.hop_size=PIL_read_mat(fid(1));
            count=0;
            for p=1:max(flq.hop_size(:,1))
                for n=1:max(flq.hop_size(:,2))
                    count=count+1;
                    flq.hop_mat(p)(n)=zeros(flq.hop_size(count,3)...
                    ,flq.hop_size(count,4),flq.hop_size(count,5));
                    for m=1:flq.hop_size(count,5)
                        flq.hop_mat(p)(n)(:,:,m)...
                        =PIL_sparse(PIL_read_mat(fid(1)),'sparse','all');
                    end
                end
            end
        case 'wan'
            wan.lat_vec=PIL_read_mat(fid(1));
            wan.rec_vec=PIL_read_mat(fid(1));
            wan.atom_type=PIL_read_mat(fid(1));
            wan.sub_lat=PIL_read_mat(fid(1));
            wan.wf_spatial=PIL_read_mat(fid(1));
            wan.wf_info=PIL_read_mat(fid(1));
            wan.uc_index=PIL_read_mat(fid(1))
            tot_wf=length(wan.wf_info(:,1));
            tot_uc=length(wan.uc_index(:,1));
            wan.Hr_mat=zeros(tot_wf,tot_wf,tot_uc);
            for n=1:tot_uc
                wan.Hr_mat(:,:,n)=PIL_read_mat(fid(1));
            end
        case 'imp'
            imp.R0_index=PIL_read_mat(fid(1));
            imp.pot_orig=PIL_read_mat(fid(1));
            imp.pot_corr=PIL_read_mat(fid(1));
        case 'hdr'
            hdr.pc_vec=PIL_read_mat(fid(1));
            hdr.pc_sublat=PIL_read_mat(fid(1));
            hdr.cc_vec=PIL_read_mat(fid(1));
            hdr.cc_sublat=PIL_read_mat(fid(1));
            hdr.slab_vec=PIL_read_mat(fid(1));
            hdr.slab_sublat=PIL_read_mat(fid(1));
            hdr.slab_atom_type=PIL_read_mat(fid(1));
            hdr.slab_pc_index=PIL_read_mat(fid(1));
            hdr.state_info=PIL_read_mat(fid(1));
            hdr.coup_site=PIL_read_mat(fid(1));
            hdr.coup_slab=PIL_read_mat(fid(1));

            tot_coup_slab=length(hdr.coup_slab(:,1));
            tot_wf=length(hdr.state_info(:,1));
            hdr.Hr_mat=zeros(tot_wf,tot_wf,tot_coup_slab);
            for n=1:tot_coup_slab
                hdr.Hr_mat(:,:,n)=PIL_read_mat(fid(1));
            end
        case 'spl'
            spl.pc_vec=PIL_read_mat(fid(1));
            spl.pc_sublat=PIL_read_mat(fid(1));
            spl.sc_vec=PIL_read_mat(fid(1));
            spl.sc_sublat=PIL_read_mat(fid(1));
            spl.sc_atom=PIL_read_mat(fid(1));
            spl.sc_pc_index=PIL_read_mat(fid(1));
            spl.state_info=PIL_read_mat(fid(1));
            spl.coup_site=PIL_read_mat(fid(1));
            spl.coup_sc=PIL_read_mat(fid(1));

            tot_coup_sc=length(spl.coup_sc(:,1));
            tot_wf=length(spl.state_info(:,1));
            spl.Hr_mat=zeros(tot_wf,tot_wf,tot_coup_sc);
            for n=1:tot_coup_sc
                spl.Hr_mat(:,:,n)=PIL_read_mat(fid(1));
            end
            // hamiltonian =============================================
        case 'ham'
            ham.lat_vec=PIL_read_mat(fid(1));
            ham.rec_vec=PIL_read_mat(fid(1));
            ham.atom_type=PIL_read_mat(fid(1));
            ham.sub_lat=PIL_read_mat(fid(1));
            ham.state_type=PIL_read_mat(fid(1));
            ham.state_info=PIL_read_mat(fid(1));
            ham.uc_index=PIL_read_mat(fid(1));

            tot_uc_index=length(ham.uc_index(:,1));
            tot_wf=length(ham.state_info(:,1));
            ham.Hr_mat=zeros(tot_wf,tot_wf,tot_uc_index);
            for n=1:tot_uc_index
                ham.Hr_mat(:,:,n)=PIL_read_mat(fid(1));
            end
            // postprocess =============================================
        case 'ban'
            ban.k_path_div=PIL_read_mat(fid(1));
            ban.k_point=PIL_read_mat(fid(1));
            ban.k_band=PIL_read_mat(fid(1));

            tot_state=length(ban.k_band(:,1));
            tot_k=length(ban.k_point(:,1));
            ban.H_k=zeros(tot_state,tot_state,tot_k);
            for n=1:tot_k
                ban.H_k(:,:,n)=PIL_read_mat(fid(1));
            end
            ban.k_vec=zeros(tot_state,tot_state,tot_k);
            for n=1:tot_k
                ban.k_vec(:,:,n)=PIL_read_mat(fid(1));
            end
        case 'ban_den'
            ban_den.site_den=PIL_read_mat(fid(1));
            ban_den.site_up_den=PIL_read_mat(fid(1));
            ban_den.site_dn_den=PIL_read_mat(fid(1));
            ban_den.orb_den=PIL_read_mat(fid(1));
            ban_den.orb_up_den=PIL_read_mat(fid(1));
            ban_den.orb_dn_den=PIL_read_mat(fid(1));
            ban_den.suborb_den=PIL_read_mat(fid(1));
        case 'dsa'
            dsa.E_dsa=PIL_read_mat(fid(1));
        case 'chn'
            chn.Chern_val=PIL_read_mat(fid(1));
            chn.V_pert=PIL_read_mat(fid(1));
            chn.k_point=PIL_read_mat(fid(1));
            chn.Fk_field=PIL_read_mat(fid(1));
        case 'zti'
            zti.Z2_val=PIL_read_mat(fid(1));
            select length(zti.Z2_val)
            case 1 //2D
                zti.V_pert=PIL_read_mat(fid(1));
                zti.n_field=PIL_read_mat(fid(1));
            case 4 //3D
                zti.Z2_plane=PIL_read_mat(fid(1));
                zti.V_pert=PIL_read_mat(fid(1));
                n_field_1=PIL_read_mat(fid(1));
                n_field_size=size(n_field_1);
                zti.n_field=zeros(n_field_size(1),n_field_size(2),3,2);
                zti.n_field(:,:,1,1)=n_field_1;
                for n=1:3
                    for m=1:2
                        if n~=1 & m~=1
                            zti.n_field(:,:,n,m)=PIL_read_mat(fid(1));
                        end
                    end
                end
            end
        end
    end
    mclose(fid(1));
    // save to .sod file #############################################
    if trim_task=='trim' then
        movefile(project_name+'_'+job_type...
        +'.tmp',project_name+'_'+job_type+'.plb');
    end
    mdelete(project_name+'_'+job_type+'.tmp');
    save(project_name+'_'+job_type+'.sod',job_type);
endfunction
