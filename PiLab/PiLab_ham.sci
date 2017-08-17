// **** Purpose ****
// This function generates ham file for postprocess calculations in PiLab
// **** Variables ****
//==== << PiLab inputs >> ====
//[ham.builder]: n x 1, string
//<= tasks that used to build Hr. Currently, PiLab support the follwing
//tasks to build Hr.
//    ['lat','hop','scc']: conventional Slaster-Koster construction
//    ['lat','hop','flq']: Slaster-Koster Floquet construction
//    ['wan']: Wannier 90 construction
//    ['wan','hdr']: Wannier 90 + dimension reduction construction
//    * the order of the tasks doesn't matter.   
//[ham.Verbosity]: 1x1, string, 'less' / 'more'
//<= if 'less', Hr related information will not be printed. 
//
//==== << PiLab outputs >> ====
//[ham.lat_vec]: 3x3, real
//=> lattice row vectros
//[ham.rec_vec]: 3x3, real
//=> reciprocal lattice row vectors
//[ham.atom_type]: tot_sublat x 1, string
//=> atom type of sublattice
//[ham.sub_lat]: tot_sublat x 3, real
//=> sublattice row vectors
//[ham.state_type]: 2x1, string
//=> first row tells the types of the states. the second
//row tells you how to understand each column in state_info
//[ham.state_info]: real, size depends on ham.state_type
//=> state info of each state, meaning of each column should
//reference to ham.state_type.
//[ham.TR_check]: 1x1, int, 0/1
//=> if 1, TRS exists. if 0, TRS not exists. 
//[ham.uc_index]: n x 3, integer
//=> unit cell index of the <n(0)|H|m(R)>
//[ham.Hr_mat]: tot_uc_index x tot_slab_wf x tot_slab_wf
//=> Hr of different R
//
// **** Version ****
// 02/15/2016 first built
// **** Comment ****

function PiLab_ham(project_name)
    disp('{ham}: starting calculation ');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{ham}: loading variables ');
    PiLab_loader(project_name,'ham','user','trim');
    load(project_name+'_ham.sod');
    disp('  => all data loaded')

    // check variables  ================================================
    disp('{ham}: checking variables ')  
    check_var=(ham.Verbosity=='more') | (ham.Verbosity=='less') 
    if check_var~=%t then
        disp('Error: PiLab_ham, ham.Verbosity must '...
        +'be ''more'' or ''less''!');
        abort;
    end
    disp('  => all variables passed')

    // core part =======================================================
    // make building tasks
    disp('{ham}: running core part ');

    task_name=['lat','hop','scc','flq','wan','hdr','spl','imp']
    builder_len=length(length(ham.builder));
    builder_id=zeros(ham.builder)
    for n=1:builder_len
        builder_id(n)=find(task_name==ham.builder(n))
    end
    builder_id=gsort(builder_id,'g','i')
    builder_text=''
    for n=1:builder_len
        builder_text=builder_text+string(builder_id(n))+'/'
    end
    disp('{ham}: generate builder IDs ');
    printf('\n   => '+strcat(repmat(' %s /',1,builder_len))+'\n'..
    ,task_name(builder_id));


    // construct Hamiltonian
    disp('{ham}: collect information from builders ');
    select builder_text
    case '1/2/3/' // conventional Slaster-Koster (lat+hop+scc)
        load(project_name+'_lat.sod');
        load(project_name+'_hop.sod');
        load(project_name+'_scc.sod');
        ham.lat_vec=lat.LatVec;
        ham.rec_vec=lat.rec_vec
        ham.sub_lat=lat.Sublat;
        ham.atom_type=lat.AtomType';
        ham.state_type(1)='Slaster-Koster'
        ham.state_type(2)='[site, identifier, l, SubOrb]'
        ham.state_info=hop.state_info(:,2:$);
        H_local=hop.onsite_E+hop.LS_mat+scc.U_mat
        [ham.uc_index,ham.Hr_mat]=..
        PIL_H_R(lat.surr_site,hop.state_info,H_local,hop.hop_mat);

        tot_state=length(ham.state_info(:,1));
        tot_uc=length(ham.uc_index(:,1));
        ham.Hr_mat=sparse(matrix(ham.Hr_mat,tot_state,tot_state*tot_uc));

    case '1/2/4/' // Slaster-Koster Floquet (lat+hop+flq)
        load(project_name+'_lat.sod');
        load(project_name+'_hop.sod');
        load(project_name+'_flq.sod');
        ham.lat_vec=lat.LatVec;
        ham.rec_vec=lat.rec_vec
        ham.sub_lat=lat.Sublat;
        ham.atom_type=lat.AtomType';
        ham.state_type(1)='Slaster-Koster Floquet'
        ham.state_type(2)='[Flq Order, site, identifier, l, SubOrb]'
        ham.state_info=flq.state_info(:,2:$);

        // generate Hr for all order
        tot_hop_state=length(hop.state_info(:,1));
        tot_flq_state=(2*flq.Order+1)*tot_hop_state;

        Hr_sub=list();
        for p=1:flq.Order+1
            [ham.uc_index,Hr_tmp]=..
            PIL_H_R(lat.surr_site,hop.state_info,flq.H_local(p),flq.hop_mat(p))
            Hr_sub(p)=Hr_tmp;
        end

        // construct flq H_R by linking different order
        tot_uc=length(ham.uc_index(:,1));
        ham.Hr_mat=zeros(tot_flq_state,tot_flq_state,tot_uc);
        // calculate upper triangle of ham.Hr_mat
        for n=1:tot_uc
            // linking different order   
            for p=1:2*flq.Order+1
                for q=p:p+flq.Order
                    r_range=(p-1)*tot_hop_state+1:p*tot_hop_state;
                    c_range=(q-1)*tot_hop_state+1:q*tot_hop_state;
                    if q <=2*flq.Order+1 then
                        if q==p & n==1 then
                            ham.Hr_mat(r_range,c_range,n)=Hr_sub(q-p+1)(:,:,n)..
                            -(p-1-flq.Order)*flq.Frequency..
                            *eye(tot_hop_state,tot_hop_state);
                        else
                            ham.Hr_mat(r_range,c_range,n)=Hr_sub(q-p+1)(:,:,n);
                        end
                    end
                end
            end       
        end 

        // calculate lower triangle of ham.Hr_mat using H(R)'=H(-R)
        for n=1:tot_uc
            mR_index=PIL_row_find(ham.uc_index,-ham.uc_index(n,:));

            ham.Hr_mat(:,:,n)=..
            triu(squeeze(ham.Hr_mat(:,:,n)))..
            +tril(squeeze(ham.Hr_mat(:,:,mR_index))')..
            -diag(diag(squeeze(ham.Hr_mat(:,:,n))));
        end

        tot_state=length(ham.state_info(:,1));
        tot_uc=length(ham.uc_index(:,1));
        ham.Hr_mat=sparse(matrix(ham.Hr_mat,tot_state,tot_state*tot_uc));
    case '5/'     // Wannier 90 (wan)
        load(project_name+'_wan.sod');
        ham.lat_vec=wan.lat_vec;
        ham.rec_vec=PIL_recip_vec(wan.lat_vec);
        ham.sub_lat=wan.sub_lat;
        ham.atom_type=wan.atom_type;
        ham.state_type(1)='Wannier-Function'
        ham.state_type(2)='[x,y,z,TR_pair]'
        ham.state_info=wan.state_info;
        ham.uc_index=wan.uc_index(:,1:3)

        ham.Hr_mat=wan.Hr_mat;        
        tot_state=length(ham.state_info(:,1));
        tot_uc=length(ham.uc_index(:,1));
        for n=1:tot_uc
            ham.Hr_mat(:,(n-1)*tot_state+1:n*tot_state)..
            =ham.Hr_mat(:,(n-1)*tot_state+1:n*tot_state)/wan.uc_index(n,4);
        end
    case '5/6/'   // Wannier 90 + dimension reduction (wan+hdr)
        load(project_name+'_wan.sod');
        load(project_name+'_hdr.sod');
        ham.lat_vec=hdr.slab_vec;
        ham.rec_vec=PIL_recip_vec(hdr.slab_vec);
        ham.sub_lat=hdr.slab_sublat(:,2:4);
        ham.atom_type=hdr.slab_atom
        ham.state_type(1)='Dimension Recuded Wannier-Function'
        ham.state_type(2)='[WF label in pc, pc index, TR_pair, (x,y,z), a3 cc index]'

        tot_state=length(hdr.state_info(:,1));
        ham.state_info=zeros(tot_state,9);
        for n=1:tot_state
            ham.state_info(n,1:4)=hdr.state_info(n,:);
            if wan.state_info(1,4)~=0 then
                ham.state_info(n,5)=PIL_row_find(hdr.state_info,..
                [wan.state_info(hdr.state_info(n,1),4),hdr.state_info(n,2:4)])
            end
            ham.state_info(n,6:8)=hdr.state_info(n,2:4)*wan.lat_vec..
            +wan.state_info(hdr.state_info(n,1),1:3);
            tmp=PIL_linexpan(ham.state_info(n,2:4)*hdr.pc_vec,hdr.cc_vec')';
            ham.state_info(n,9)=floor(tmp(3)+1e-4);
        end
        ham.uc_index=hdr.uc_index;
        ham.Hr_mat=hdr.Hr_mat;
        clear hdr wan
    case '5/7/' // wannier90 + super lattice (wan+spl)
        load(project_name+'_wan.sod');
        load(project_name+'_spl.sod');
        ham.lat_vec=spl.sc_vec;
        ham.rec_vec=PIL_recip_vec(spl.sc_vec);
        ham.sub_lat=spl.sc_sublat(:,2:4);
        ham.atom_type=spl.sc_atom
        ham.state_type(1)='Superlattice Wannier-Function'
        ham.state_type(2)='[WF label in pc, pc index, TR_pair, (x,y,z) ]'

        tot_state=length(spl.state_info(:,1));
        ham.state_info=zeros(tot_state,8);
        for n=1:tot_state
            ham.state_info(n,1:4)=spl.state_info(n,:);
                if wan.state_info(1,4)~=0 then
                ham.state_info(n,5)=PIL_row_find(spl.state_info,..
                [wan.state_info(spl.state_info(n,1),4),spl.state_info(n,2:4)])
            end
            ham.state_info(n,5:7)=spl.state_info(n,2:4)*wan.lat_vec..
            +wan.state_info(spl.state_info(n,1),1:3);
        end
        ham.uc_index=spl.uc_index;
        ham.Hr_mat=spl.Hr_mat;
    case '5/8/'  // wan --> imp
        load(project_name+'_wan.sod');
        ham.lat_vec=wan.lat_vec;
        ham.rec_vec=PIL_recip_vec(wan.lat_vec);
        ham.sub_lat=wan.sub_lat;
        ham.atom_type=wan.atom_type;
        ham.state_type(1)='Wannier-Function'
        ham.state_type(2)='[sublat, norm_err, TR_pair]'
        ham.state_info=wan.state_info(:,[1,5,6])
        ham.uc_index=wan.uc_index(:,1:3)

        ham.Hr_mat=wan.Hr_mat;        
        tot_state=length(ham.state_info(:,1));
        tot_uc=length(ham.uc_index(:,1));
        for n=1:tot_uc
            ham.Hr_mat(:,(n-1)*tot_state+1:n*tot_state)..
            =ham.Hr_mat(:,(n-1)*tot_state+1:n*tot_state)/wan.uc_index(n,4);
        end

        load(project_name+'_imp.sod');
        if imp.HrRead~='wan'
            disp('Error: PiLab_ham, ham.builder is inconsistent'..
            +' with imp.HrRead !');
            abort
        end
        ham.Hr_mat(:,(imp.R0_index-1)*tot_state+1:imp.R0_index*tot_state)=..
        ham.Hr_mat(:,(imp.R0_index-1)*tot_state+1:imp.R0_index*tot_state)..
        +sparse(diag(imp.pot_corr))
    case '5/6/8/' // wan -> hdr -> imp
        load(project_name+'_wan.sod');
        load(project_name+'_hdr.sod');
        ham.lat_vec=hdr.slab_vec;
        ham.rec_vec=PIL_recip_vec(hdr.slab_vec);
        ham.sub_lat=hdr.slab_sublat(:,2:4);
        ham.atom_type=hdr.slab_atom
        ham.state_type(1)='Dimension Recuded Wannier-Function'
        ham.state_type(2)='[WF label in pc, pc index, TR_pair, (x,y,z), cc a3 index]'

        tot_state=length(hdr.state_info(:,1));
        ham.state_info=zeros(tot_state,9);
        for n=1:tot_state
            ham.state_info(n,1:4)=hdr.state_info(n,:);
            if wan.state_info(1,4)~=0 then
                ham.state_info(n,5)=PIL_row_find(hdr.state_info,..
                [wan.state_info(hdr.state_info(n,1),4),hdr.state_info(n,2:4)])
            end
            ham.state_info(n,6:8)=hdr.state_info(n,2:4)*wan.lat_vec..
            +wan.state_info(hdr.state_info(n,1),1:3);
            tmp=PIL_linexpan(ham.state_info(n,2:4)*hdr.pc_vec,hdr.cc_vec')';
            ham.state_info(n,9)=floor(tmp(3)+1e-4);
        end
        ham.uc_index=hdr.uc_index;
        ham.Hr_mat=hdr.Hr_mat;

        load(project_name+'_imp.sod');
        if imp.HrRead~='hdr'
            disp('Error: PiLab_ham, ham.builder is inconsistent'..
            +' with imp.HrRead !');
            abort
        end
        tot_state=length(ham.state_info(:,1));
        ham.Hr_mat(:,(imp.R0_index-1)*tot_state+1:imp.R0_index*tot_state)=..
        ham.Hr_mat(:,(imp.R0_index-1)*tot_state+1:imp.R0_index*tot_state)..
        +sparse(diag(imp.pot_corr));
    case '5/7/8/' // wan --> spl --> imp
        load(project_name+'_wan.sod');
        load(project_name+'_spl.sod');
        ham.lat_vec=spl.sc_vec;
        ham.rec_vec=PIL_recip_vec(spl.sc_vec);
        ham.sub_lat=spl.sc_sublat(:,2:4);
        ham.atom_type=spl.sc_atom
        ham.state_type(1)='Superlattice Wannier-Function'
        ham.state_type(2)='[WF label in pc, pc index, TR_pair, (x,y,z) ]'

        tot_state=length(spl.state_info(:,1));
        ham.state_info=zeros(tot_state,8);
        for n=1:tot_state
            ham.state_info(n,1:4)=spl.state_info(n,:);
                if wan.state_info(1,4)~=0 then
                ham.state_info(n,5)=PIL_row_find(spl.state_info,..
                [wan.state_info(spl.state_info(n,1),4),spl.state_info(n,2:4)])
            end
            ham.state_info(n,5:7)=spl.state_info(n,2:4)*wan.lat_vec..
            +wan.state_info(spl.state_info(n,1),1:3);
        end
        ham.uc_index=spl.uc_index;
        ham.Hr_mat=spl.Hr_mat;
        
        load(project_name+'_imp.sod');
        if imp.HrRead~='spl'
            disp('Error: PiLab_ham, ham.builder is inconsistent'..
            +' with imp.HrRead !');
            abort
        end
        tot_state=length(ham.state_info(:,1));
        ham.Hr_mat(:,(imp.R0_index-1)*tot_state+1:imp.R0_index*tot_state)=..
        ham.Hr_mat(:,(imp.R0_index-1)*tot_state+1:imp.R0_index*tot_state)..
        +sparse(diag(imp.pot_corr))
    else
        disp('Error: PiLab_ham, ham.builder are not correct builders!')
        abort
    end
    disp('  => all information collected')

    // check hermitian 
    disp('{ham}: checking hermitianity' )
    [Hk,err]=PIL_Hk_R([1,2,3],ham.lat_vec,ham.uc_index,[],ham.Hr_mat);
    if err==1 then
        abort
    else
        disp('  => passed !');
    end
    // output information ==============================================
    disp('{ham}: output information ')
    fid=mopen(project_name+'_ham.plb','a+'); 
    PIL_print_mat('ham.lat_vec, @f:f, lattice row vectors',..
    ham.lat_vec,'r',fid(1));
    PIL_print_mat('ham.rec_vec, @f:f, reciprocal lattice row vectors',..
    ham.rec_vec,'r',fid(1));
    PIL_print_mat('ham.atom_type, @f:f, atom type of sublattice',..
    ham.atom_type,'s',fid(1));
    PIL_print_mat('ham.sub_lat, @f:f, sublattice row vectors',..
    ham.sub_lat,'r',fid(1));
    PIL_print_mat('ham.state_type, @f:f, type of states',..
    ham.state_type,'s',fid(1));
    PIL_print_mat('ham.state_info, @f:f, information of each state',..
    ham.state_info,'r',fid(1),'off');
    PIL_print_mat('ham.uc_index, @f:f, unit cell index in H(R)',..
    ham.uc_index,'i',fid(1));
    mclose(fid)
    if ham.Verbosity=='more' then
        tot_uc=length(ham.uc_index(:,1))
        tot_state=length(ham.state_info(:,1))
        printf('\n   => writing real part of Hr to '+..
        project_name+'_ham_Hr_real.dat\n\n');
        fprintfMat(project_name+'_ham_Hr_real.dat',..
        (real(full(ham.Hr_mat)))'..
        ,'%f','real part of Hr, size=['..
        +string(tot_uc*tot_state)+','+string(tot_state)+']');
        
        printf('   => writing imag part of Hr to '+..
        project_name+'_ham_Hr_imag.dat\n');
        fprintfMat(project_name+'_ham_Hr_imag.dat',..
        (imag(full(ham.Hr_mat)))'..
        ,'%f','imag part of Hr, size=['..
        +string(tot_uc*tot_state)+','+string(tot_state)+']');
    end
    
    // finishing program ===============================================
    save(project_name+'_ham.sod','ham');
    disp('{ham}: finishing calculation ');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction
