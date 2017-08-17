// **** Purpose ****
// PiLab hamiltonian dimension reduction function
// **** Variables ****
//==== << PiLab inputs >> ====
//[spl.SuperCell]: 3x3, real
//<= the conventional cell row vecvtors, in unit of the WF lattice vectors 
//[spl.Task]: 1x1, string, 'str' / 'run'
//<= 'str': output structure file only. 'run': run whole module.  
//[spl.Verbosity]: 1x1, string, 'less' / 'more'
//<= if 'less', H_r related information will not br printed.
//
//==== << PiLab outputs >> ====
//[spl.pc_vec]: 3x3, real
//=> lattice row vectors of primitive cell (w90 input)
//[spl.pc_sublat]: nx3, real
//=> sublattice row vectors of primitive sublattice (w90 input) 
//[spl.pc_atom]: nx1, string
//=> atom type of primitive cell
//
//[spl.sc_vec]: 3x3, real
//=> lattice row vectors of sc unitcell
//[spl.sc_sublat]: tot_sc_sublat x 5, real
//=> sublat index in pc(1), x,y,z, proj on a3
//[spl.sc_atom]: tot_sc_sublat x 1, string
//=> atom type of sc sublattice
//
//[spl.sc_unit]: nx10, real
//=> primitive unit cell inside the sc
//[(b=1,n1,n2,n3),x,y,z,expan in sc]
//
//[spl.coup_unit]: nx7, int  
//=> pc coupled to the R=0 sc
//[(n1,n2,n3 in pc),(n1, n2, n3 in sc), uc_deg]
//
//[spl.state_info]: tot_sc_wf x 4, integer 
//=> state info of sc WFs
//[(wf label in pc), (n1, n2, n3 in pc)]
//
//[spl.uc_index]: nx3, integer
//=> R of Hr_mat index in sc lattice vector
//
//[spl.Hr_mat]:  tot_sc_WF x tot_sc_WF x tot_coup_sc
//=> Hr of the sc structure, sparse format
//
// **** Version ****
// 02/10/2016 first built
// **** Comment ****

function PiLab_spl(project_name)
    disp('{spl}: starting calculation ...');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{spl}: loading variables ...');
    PiLab_loader(project_name,'spl','user','trim');
    load(project_name+'_spl.sod');
    load(project_name+'_wan.sod');
    disp('  => all data loaded')

    // check variables  ================================================
    disp('{spl}: checking variables ...')  
    check_var=(length(spl.SuperCell(:,1)==3)) & (length(spl.SuperCell(1,:)==3));
    if check_var~=%t then
        disp('Error: PiLab_spl, spl.SuperCell must '...
        +'be 3x3 matrix !');
        abort;
    end
    check_var=(sum(spl.SuperCell-round(spl.SuperCell))==0)
    if check_var~=%t then
        disp('Error: PiLab_spl, spl.SuperCell must '...
        +'be a 3x3 integer matrix !');
        abort;
    end
    check_var=(spl.Task=='plot') | (spl.Task=='run') 
    if check_var~=%t then
        disp('Error: PiLab_spl, spl.Task must '...
        +'be ''plot'' or ''run''!');
        abort;
    end
    check_var=(spl.Verbosity=='more') | (spl.Verbosity=='less') 
    if check_var~=%t then
        disp('Error: PiLab_spl, spl.Verbosity must '...
        +'be ''more'' or ''less''!');
        abort;
    end
    disp('  => all variables passed')

    // core part =======================================================
    disp('{spl}: running core part ...');

    // eliminate unwanted sites  ---------------------------------------
    disp('  => eliminating unwanted sites');
    spl.pc_vec=wan.lat_vec;
    spl.pc_sublat=wan.sub_lat;
    spl.pc_atom=wan.atom_type;
    tot_pc_sub=length(spl.pc_sublat(:,1))
    tot_pc_wf=length(wan.wf_spatial(:,1));

    // generate conventional cell along assigned direction -------------
    disp('  => building spuer cell')
    [spl.sc_sublat]=PIL_conv_cell_vec(spl.pc_vec,spl.pc_sublat,spl.SuperCell);
    spl.sc_sublat=clean(spl.sc_sublat,1e-5);
    spl.sc_vec=clean(PIL_red_cart_conv(spl.pc_vec,spl.SuperCell,'red'),1e-5);
    spl.sc_atom=spl.pc_atom(spl.sc_sublat(:,1));

    // output xsf file for visualization -------------------------------
    disp('  => output structure files')

    PIL_crystal_xsf(project_name+'_spl_pc',spl.pc_vec,..
    spl.pc_atom,spl.pc_sublat(:,1:3))

    PIL_crystal_xsf(project_name+'_spl_cc',spl.sc_vec,..
    spl.sc_atom,spl.sc_sublat(:,5:7))

    select spl.Task
    case 'plot'
        disp('      * Warning:') 
        disp('        structure ouput, calculation not completed!');
        save(project_name+'_spl.sod','spl');
        fid=mopen(project_name+'_spl.plb','a+'); 
        PIL_print_mat('spl.pc_vec, @f:f, lattice vectors of primitive cell',..
        spl.pc_vec,'r',fid(1));
        PIL_print_mat('spl.pc_sublat, @f:f, sublattices of primitive cell',..
        spl.pc_sublat,'r',fid(1));
        PIL_print_mat('spl.pc_atom, @f:f, atom types of primitive cell',..
        spl.pc_atom,'s',fid(1));

        PIL_print_mat('spl.sc_vec, @f:f, lattice vectors of conventional cell',..
        spl.sc_vec,'r',fid(1));
        PIL_print_mat('spl.sc_sublat, @f:f, sublattice of conventional cell'..
        +' [index in pc(4),x, y, z, expan in cc(3)]',..
        spl.sc_sublat,'r',fid(1),'off');
        PIL_print_mat('spl.sc_atom, @f:f, atom types of super cell',..
        spl.sc_atom,'s',fid(1));
        mclose(fid);
        abort
    end

    // generate sc wf locations --------------------------------------
    disp('  => searching for coupled unit cell');
    // generate pc unit in cc and sc
    // imagine there are tot_state in a single pc unit cell
    // spl.sc_unit: [b(always 1), n1,n2,n3 in pc, x,y,z, sc expan]
    [spl.sc_unit]=clean(PIL_conv_cell_vec(spl.pc_vec,[0,0,0],spl.SuperCell),1e-5);

    // search for all coupled units ------------------------------------
    // coup_unit_tmp:[n1,n2,n3,uc_deg]
    tot_pc_uc_index=length(wan.uc_index(:,1));
    tot_sc_unit=length(spl.sc_unit(:,1));
    coup_unit_tmp=zeros(tot_pc_uc_index*tot_sc_unit,4);
    for n=1:tot_sc_unit
        coup_unit_tmp(tot_pc_uc_index*(n-1)+1:tot_pc_uc_index*n,:)=..
        wan.uc_index+cat(2,repmat(spl.sc_unit(n,2:4),tot_pc_uc_index,1)..
        ,zeros(tot_pc_uc_index,1));
    end
    coup_unit_tmp=gsort(coup_unit_tmp,'lr','i');

    // eliminate repeated or out-of-range coup_unit
    // spl.coup_unit:[pc index x3, sc index x3, uc_deg]
    spl.coup_unit=zeros(length(coup_unit_tmp(:,1)),7);
    count=0;
    for n=1:length(coup_unit_tmp(:,1))
        if n~=1 & sum(abs(coup_unit_tmp(n,:)-coup_unit_tmp(n-1,:))) <=1e-4 then
            continue;
        end
        sc_index=PIL_linexpan(coup_unit_tmp(n,1:3)*spl.pc_vec,spl.sc_vec');
        sc_index=(sc_index+1e-4)'; // numerical error to integer
        count=count+1;
        spl.coup_unit(count,1:3)=coup_unit_tmp(n,1:3);
        spl.coup_unit(count,4:6)=floor(sc_index);
        spl.coup_unit(count,7)=coup_unit_tmp(n,4);

    end
    spl.coup_unit=spl.coup_unit(1:count,:);
    
    // check result
    tmp1=spl.coup_unit(PIL_row_find(spl.coup_unit(:,4:6),[0,0,0]),1:3);
    if PIL_equal(gsort(tmp1,'lr','i'),gsort(spl.sc_unit(:,2:4),'lr','i')) then

    else
        disp('Error: PiLab_spl, spl.coup_unit doesn''t match spl.sc_unit!');
        abort
    end
    clear coup_unit_tmp;

    // generate coupled sc index list --------------------------------
    tot_coup_unit=length(spl.coup_unit(:,1));
    spl.uc_index=zeros(tot_coup_unit,3);
    coup_sc_tmp=gsort(spl.coup_unit(:,4:6),'lr','i');
    spl.uc_index(1,:)=coup_sc_tmp(1,:);
    count=1;
    for n=2:tot_coup_unit
        if sum(abs(coup_sc_tmp(n,:)-coup_sc_tmp(n-1,:)))<1e-4 then
            continue;
        end
        count=count+1;
        spl.uc_index(count,:)=coup_sc_tmp(n,:);
    end
    spl.uc_index=spl.uc_index(1:count,:);

    // build Hamiltonian -----------------------------------------------
    disp('  => building Hamiltonian');
    // sc_wf_list:[wf_label, label of sc_unit]
    tot_sc_wf=tot_pc_wf*tot_sc_unit;
    sc_wf_list=matrix(1:tot_sc_wf,tot_pc_wf,tot_sc_unit);
    spl.Hr_mat=spzeros(tot_sc_wf,tot_sc_wf*length(spl.uc_index(:,1)));
    spl.state_info=zeros(tot_sc_wf,4);
    for n=1:tot_sc_unit
        n_wf=sc_wf_list(:,n);
        spl.state_info((n-1)*tot_pc_wf+1:n*tot_pc_wf,:)=..
        cat(2,[1:tot_pc_wf]',repmat(spl.sc_unit(n,2:4),tot_pc_wf,1));
        for m=1:tot_coup_unit
            // locate all index
            R_sc_idx=PIL_row_find(spl.uc_index,spl.coup_unit(m,4:6));

            m_idx=PIL_row_find(spl.sc_unit(:,8:10),..
            PIL_linexpan(spl.coup_unit(m,1:3)*spl.pc_vec..
            -spl.coup_unit(m,4:6)*spl.sc_vec,spl.sc_vec')');

            R_pc_idx=PIL_row_find(wan.uc_index(:,1:3),..
            spl.coup_unit(m,1:3)-spl.sc_unit(n,2:4));

            if R_pc_idx~=[] & m_idx~=[] & R_sc_idx~=[] then
                m_wf=sc_wf_list(:,m_idx);
                spl.Hr_mat(n_wf,m_wf+(R_sc_idx-1)*tot_sc_wf)=..
                wan.Hr_mat(:,(R_pc_idx-1)*tot_pc_wf+1:(R_pc_idx)*tot_pc_wf)..
                /wan.uc_index(R_pc_idx,4);     
            end
        end
    end

    // check if coupled scs are enough for hermitian -----------------
    disp('  => checking Hr Hermiticity')
    [Hk,err]=PIL_Hk_R([1,2,3],spl.sc_vec,spl.uc_index,[],spl.Hr_mat);
    if err==1 then
        disp('Error: PiLab_spl, Hr is not hermitian !')
        abort
    else
        disp('      passed!')
    end

    // output information ==============================================
    disp('{spl}: output information ...')
    fid=mopen(project_name+'_spl.plb','a+'); 

    PIL_print_mat('spl.pc_vec, @f:f, lattice vectors of primitive cell',..
    spl.pc_vec,'r',fid(1));
    PIL_print_mat('spl.pc_sublat, @f:f, sublattices of primitive cell',..
    spl.pc_sublat,'r',fid(1));
    PIL_print_mat('spl.pc_atom, @f:f, atom type of primitive cell',..
    spl.pc_atom,'s',fid(1));

    PIL_print_mat('spl.sc_vec, @f:f, lattice vectors of super cell',..
    spl.sc_vec,'r',fid(1));
    PIL_print_mat('spl.sc_sublat, @f:f, sublattice of super cell'..
    +' [index in pc(4),x, y, z, expan in cc(3)]',..
    spl.sc_sublat,'r',fid(1),'off');
    PIL_print_mat('spl.sc_atom, @f:f, atom type of super cell',..
    spl.sc_atom,'s',fid(1));

    PIL_print_mat('spl.sc_unit, @f:f, pc in sc,'..
    +'[b=1,n1,n2,n3,x,y,z,expan in sc]',spl.sc_unit,'r',fid(1));

    PIL_print_mat('spl.coup_unit, @f:f, pc coupled to the R=0 sc,'..
    +'[(n1,n2,n3 in pc),(n1, n2, n3 in sc), uc_deg]',spl.coup_unit,'i',fid(1));

    PIL_print_mat('spl.state_info, @f:f, state info of sc WFs,'..
    +'[(wf label in pc), (n1, n2, n3 in pc)]',spl.state_info,'i',fid(1));

    if spl.Verbosity=='more' then
        PIL_print_mat('spl.uc_index, @f:f, sc cells coupled to R=0 sc cell'..
        +' [n1,n2,n3 in sc]',spl.uc_index,'i',fid(1));
        for n=1:length(spl.uc_index(:,1))
            PIL_print_mat('spl.Hr_mat(:,(n-1)*tot_sc_wf+1:n*tot_sc_wf)'+..
            ', @as:as, Hr matrix of '+string(n)+'-th sc_index'..
            ,spl.Hr_mat(:,(n-1)*tot_sc_wf+1:n*tot_sc_wf)..
            ,'sp',fid(1));
        end
    end
    mclose(fid(1))

    // finishing program ===============================================
    save(project_name+'_spl.sod','spl');
    disp('{spl}: finishing calculation ...');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction
