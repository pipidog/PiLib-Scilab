// **** Purpose ****
// PiLab hamiltonian dimension reduction function
// **** Variables ****
//==== << PiLab inputs >> ====
//[hdr.ConvCell]: 3x3, real
//<= the conventional cell row vecvtors, in unit of the WF lattice vectors 
//[hdr.TotCC]: 1x1, integer
//<= number of conventional cell along a3 direction
//[hdr.Task]: 1x1, string, 'str' / 'run'
//<= 'str': output structure file only. 'run': run whole module.  
//[hdr.Verbosity]: 1x1, string, 'less' / 'more'
//<= if 'less', H_r related information will not br printed.
//
//==== << PiLab outputs >> ====
//[hdr.pc_vec]: 3x3, real
//=> lattice row vectors of primitive cell (w90 input)
//[hdr.pc_sublat]: nx3, real
//=> sublattice row vectors of primitive sublattice (w90 input) 
//[hdr.pc_atom]: nx1, string
//=> atom type of primitive cell
//
//[hdr.cc_vec]: 3x3, real
//=> lattic row vector of conventional cell 
//[hdr.cc_sublat]: tot_cc_sublat x 10, real
//=> index in pc(4),x, y, z, expan in cc(3)
//[hdr.cc_atom]: tot_cc_sublat x 1, string
//=> atom type of conventional cell
//
//[hdr.slab_vec]: 3x3, real
//=> lattice row vectors of slab unitcell
//[hdr.slab_sublat]: tot_slab_sublat x 5, real
//=> sublat index in pc(1), x,y,z, proj on a3
//[hdr.slab_atom]: tot_slab_sublat x 1, string
//=> atom type of slab sublattice
//
//[hdr.cc_unit]: nx10, real
//=> primitive unit cell inside the conevntional cell
//  [(b=1,n1,n2,n3),x,y,z,expan in cc]
//[hdr.slab_unit]: nx10, real
//=> primitive unit cell inside the slab cell
//  [(b=1,n1,n2,n3),x,y,z,expan in sc]
//  
//[hdr.coup_unit]: nx7, int  
//=> pc coupled to the R=0 slab
//  [(n1,n2,n3 in pc),(n1, n2, n3 in slab), uc_deg]
//
//[hdr.state_info]: tot_slab_wf x 4, integer 
//=> state info of slab WFs
//   [(wf label in pc), (n1, n2, n3 in pc)]
//
//[hdr.uc_index]: nx3, integer
//=> R of Hr_mat index in slab lattice vector
//
//[hdr.Hr_mat]:  tot_slab_WF x tot_slab_WF x tot_coup_slab
//=> Hr of the slab structure, sparse format
//
//
// **** Version ****
// 02/10/2016 first built
// **** Comment ****

function PiLab_hdr(project_name)
    disp('{hdr}: starting calculation ');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{hdr}: loading variables ');
    PiLab_loader(project_name,'hdr','user','trim');
    load(project_name+'_hdr.sod');
    load(project_name+'_wan.sod');
    disp('  => all data loaded')

    // check variables  ================================================
    disp('{hdr}: checking variables ')  
    check_var=(length(hdr.ConvCell(:,1)==3)) & (length(hdr.ConvCell(1,:)==3));
    if check_var~=%t then
        disp('Error: PiLab_hdr, hdr.ConvCell must '...
        +'be 3x3 matrix !');
        abort;
    end
    check_var=(sum(hdr.ConvCell-round(hdr.ConvCell))==0)
    if check_var~=%t then
        disp('Error: PiLab_hdr, hdr.ConvCell must '...
        +'be a 3x3 integer matrix !');
        abort;
    end
    check_var=(hdr.Task=='plot') | (hdr.Task=='run') 
    if check_var~=%t then
        disp('Error: PiLab_hdr, hdr.Task must '...
        +'be ''plot'' or ''run''!');
        abort;
    end
    check_var=(hdr.Verbosity=='more') | (hdr.Verbosity=='less') 
    if check_var~=%t then
        disp('Error: PiLab_hdr, hdr.Verbosity must '...
        +'be ''more'' or ''less''!');
        abort;
    end
    disp('  => all variables passed')

    // core part =======================================================
    disp('{hdr}: running core part ');

    // eliminate unwanted sites  ---------------------------------------
    disp('  => eliminating unwanted sites');
    hdr.pc_vec=wan.lat_vec;
    hdr.pc_sublat=wan.sub_lat;
    hdr.pc_atom=wan.atom_type;
    tot_pc_sub=length(hdr.pc_sublat(:,1))
    tot_pc_wf=length(wan.wf_spatial(:,1));

    // generate conventional cell along assigned direction -------------
    disp('  => building conventional cell')
    [hdr.cc_sublat]=PIL_conv_cell_vec(hdr.pc_vec,hdr.pc_sublat,hdr.ConvCell);
    hdr.cc_sublat=clean(hdr.cc_sublat,1e-5);
    hdr.cc_vec=clean(PIL_red_cart_conv(hdr.pc_vec,hdr.ConvCell,'red'),1e-5);
    hdr.cc_atom=hdr.pc_atom(hdr.cc_sublat(:,1));

    // generate slab structure -----------------------------------------
    disp('  => building slab structure')
    [hdr.slab_vec,hdr.slab_sublat]..
    =PIL_slab_str(hdr.cc_vec,hdr.cc_sublat,hdr.TotCC,[],[]);
    hdr.slab_atom=wan.atom_type(hdr.slab_sublat(:,1));


    // output xsf file for visualization -------------------------------
    disp('  => output structure files')
    atom_label=list()
    atom_label(1)=hdr.pc_atom;
    atom_label(2)=string(hdr.cc_sublat(:,1));
    atom_label(3)=string(hdr.slab_sublat(:,1));
    for n=1:length(hdr.pc_sublat(:,1))
        for m=2:3
            atom_label(m)(atom_label(m)==string(n))=hdr.pc_atom(n);
        end
    end

    PIL_crystal_xsf(project_name+'_hdr_pc',hdr.pc_vec,..
    atom_label(1),hdr.pc_sublat(:,1:3))

    PIL_crystal_xsf(project_name+'_hdr_cc',hdr.cc_vec,..
    atom_label(2),hdr.cc_sublat(:,5:7))
    
    slab_vec_tmp=hdr.slab_vec;
    slab_vec_tmp(3,:)=(hdr.TotCC+1)/hdr.TotCC*hdr.slab_vec(3,:);
    PIL_crystal_xsf(project_name+'_hdr_slab',slab_vec_tmp,..
    atom_label(3),hdr.slab_sublat(:,2:4)); // use larger slab_vec for plot

    select hdr.Task
    case 'plot'
        disp('      * Warning:') 
        disp('        structure ouput, calculation not completed!');
        save(project_name+'_hdr.sod','hdr');
        fid=mopen(project_name+'_hdr.plb','a+'); 
        PIL_print_mat('hdr.pc_vec, @f:f, lattice vectors of primitive cell',..
        hdr.pc_vec,'r',fid(1));
        PIL_print_mat('hdr.pc_sublat, @f:f, sublattices of primitive cell',..
        hdr.pc_sublat,'r',fid(1));
        PIL_print_mat('hdr.pc_atom, @f:f, atom type of primitive cell'..
        ,hdr.pc_atom,'s',fid(1))  

        PIL_print_mat('hdr.cc_vec, @f:f, lattice vectors of conventional cell',..
        hdr.cc_vec,'r',fid(1));
        PIL_print_mat('hdr.cc_sublat, @f:f, sublattice of conventional cell'..
        +' [index in pc(4),x, y, z, expan in cc(3)]',..
        hdr.cc_sublat,'r',fid(1),'off');
        PIL_print_mat('hdr.cc_atom, @f:f, atom type of conventional cell'..
        ,hdr.cc_atom,'s',fid(1))  

        PIL_print_mat('hdr.slab_vec, @f:f, lattice vectors of slab cell',..
        hdr.slab_vec,'r',fid(1));
        PIL_print_mat('hdr.slab_sublat, @f:f, sublat of slab'..
        +' [sublat index in pc(1), x,y,z, proj on unit a3]',..
        hdr.slab_sublat,'r',fid(1),'off');  
        PIL_print_mat('hdr.slab_atom, @f:f, atom type of slab sublattice'..
        ,hdr.slab_atom,'s',fid(1))  
        mclose(fid);
        abort
    end

    // generate slab wf locations --------------------------------------
    disp('  => searching for coupled unit cell');
    // generate pc unit in cc and slab
    // imagine there are tot_state in a single pc unit cell
    // hdr.cc_unit: [b(always 1), n1,n2,n3 in pc, x,y,z, sc expan]
    // hdr.slab_unit: [b(always 1), n1,n2,n3 in pc, x,y,z, sc expan]
    [hdr.cc_unit]=clean(PIL_conv_cell_vec(hdr.pc_vec,[0,0,0],hdr.ConvCell),1e-5);
    hdr.slab_vec=hdr.ConvCell;
    hdr.slab_vec(3,:)=hdr.TotCC*hdr.slab_vec(3,:);
    [hdr.slab_unit]=clean(PIL_conv_cell_vec(hdr.pc_vec,[0,0,0],hdr.slab_vec),1e-5);
    hdr.slab_vec=hdr.slab_vec*hdr.pc_vec;
    
    // search for all coupled units ------------------------------------
    // coup_unit_tmp:[n1,n2,n3,uc_deg]
    tot_pc_uc_index=length(wan.uc_index(:,1));
    tot_slab_unit=length(hdr.slab_unit(:,1));
    coup_unit_tmp=zeros(tot_pc_uc_index*tot_slab_unit,4);
    for n=1:tot_slab_unit
        coup_unit_tmp(tot_pc_uc_index*(n-1)+1:tot_pc_uc_index*n,:)=..
        wan.uc_index+cat(2,repmat(hdr.slab_unit(n,2:4),tot_pc_uc_index,1)..
        ,zeros(tot_pc_uc_index,1));
    end
    coup_unit_tmp=gsort(coup_unit_tmp,'lr','i');

    // eliminate repeated or out-of-range coup_unit
    // hdr.coup_unit:[pc index x3, slab index x3, uc_deg]
    hdr.coup_unit=zeros(length(coup_unit_tmp(:,1)),7);
    count=0;
    for n=1:length(coup_unit_tmp(:,1))
        if n~=1 & sum(abs(coup_unit_tmp(n,:)-coup_unit_tmp(n-1,:))) <=1e-4 then
            continue;
        end
        slab_index=PIL_linexpan(coup_unit_tmp(n,1:3)*hdr.pc_vec,hdr.slab_vec');
        slab_index=(slab_index+1e-5)'; // numerical error to integer
        if slab_index(3) >= 0 & slab_index(3) < 1 then
            count=count+1;
            hdr.coup_unit(count,1:3)=coup_unit_tmp(n,1:3);
            hdr.coup_unit(count,4:6)=floor(slab_index);
            hdr.coup_unit(count,7)=coup_unit_tmp(n,4);
        end
    end
    hdr.coup_unit=hdr.coup_unit(1:count,:);
    
//    // check if pc in slab [0,0,0] equal to the results in hdr.slab_unit 
//    tmp1=hdr.coup_unit(PIL_row_find(hdr.coup_unit(:,4:6),[0,0,0]),1:3);
//    disp((gsort(tmp1,'lr','i')));
//    disp((gsort(hdr.slab_unit(:,2:4),'lr','i')));
//    abort
//    if PIL_equal(gsort(tmp1,'lr','i'),gsort(hdr.slab_unit(:,2:4),'lr','i')) then
//
//    else
//        disp('Error: PiLab_hdr, hdr.coup_unit doesn''t match hdr.slab_unit!');
//        abort
//    end
//    if find(hdr.coup_unit(:,6)~=0)~=[] then
//        disp('Error: PiLab_hdr, hdr.coup_unit has non-zero slab a3 component!');
//    end
//    clear coup_unit_tmp;
    
    // generate coupled slab index list --------------------------------
    tot_coup_unit=length(hdr.coup_unit(:,1));
    hdr.uc_index=zeros(tot_coup_unit,3);
    coup_slab_tmp=gsort(hdr.coup_unit(:,4:6),'lr','i');
    hdr.uc_index(1,:)=coup_slab_tmp(1,:);
    count=1;
    for n=2:tot_coup_unit
        if sum(abs(coup_slab_tmp(n,:)-coup_slab_tmp(n-1,:)))<1e-4 then
            continue;
        end
        count=count+1;
        hdr.uc_index(count,:)=coup_slab_tmp(n,:);
    end
    hdr.uc_index=hdr.uc_index(1:count,:);

    // build Hamiltonian -----------------------------------------------
    disp('  => building Hamiltonian');
    // slab_wf_list:[wf_label, label of slab_unit]
    tot_slab_wf=tot_pc_wf*tot_slab_unit;
    slab_wf_list=matrix(1:tot_slab_wf,tot_pc_wf,tot_slab_unit);
    hdr.Hr_mat=spzeros(tot_slab_wf,tot_slab_wf*length(hdr.uc_index(:,1)));
    hdr.state_info=zeros(tot_slab_wf,4);
    printf('      total primitive cell in slab=%i\n', tot_slab_unit);
    for n=1:tot_slab_unit
        printf('        constructing Hamiltonian n=%i\n',n);
        n_wf=slab_wf_list(:,n);
        hdr.state_info((n-1)*tot_pc_wf+1:n*tot_pc_wf,:)=..
        cat(2,[1:tot_pc_wf]',repmat(hdr.slab_unit(n,2:4),tot_pc_wf,1));
        for m=1:tot_coup_unit
            // locate all index
            R_slab_idx=PIL_row_find(hdr.uc_index,hdr.coup_unit(m,4:6));

            m_idx=PIL_row_find(hdr.slab_unit(:,8:10),..
            PIL_linexpan(hdr.coup_unit(m,1:3)*hdr.pc_vec..
            -hdr.coup_unit(m,4:6)*hdr.slab_vec,hdr.slab_vec')');

            R_pc_idx=PIL_row_find(wan.uc_index(:,1:3),..
            hdr.coup_unit(m,1:3)-hdr.slab_unit(n,2:4));

            if R_pc_idx~=[] & m_idx~=[] & R_slab_idx~=[] then
                m_wf=slab_wf_list(:,m_idx);
                hdr.Hr_mat(n_wf,m_wf+(R_slab_idx-1)*tot_slab_wf)=..
                wan.Hr_mat(:,(R_pc_idx-1)*tot_pc_wf+1:(R_pc_idx)*tot_pc_wf)..
                /wan.uc_index(R_pc_idx,4);     
            end
        end
    end
    // check if coupled slabs are enough for hermitian -----------------
    disp('  => checking Hr Hermiticity')
    [Hk,err]=PIL_Hk_R([1,2,3],hdr.slab_vec,hdr.uc_index,[],hdr.Hr_mat);
    if err==1 then
        disp('Error: PiLab_hdr, Hr is not hermitian !')
        abort
    else
        disp('      passed!')
    end

    // output information ==============================================
    disp('{hdr}: output information ')
    fid=mopen(project_name+'_hdr.plb','a+'); 

    PIL_print_mat('hdr.pc_vec, @f:f, lattice vectors of primitive cell',..
    hdr.pc_vec,'r',fid(1));
    PIL_print_mat('hdr.pc_sublat, @f:f, sublattices of primitive cell',..
    hdr.pc_sublat,'r',fid(1));
    PIL_print_mat('hdr.pc_atom, @f:f, atom type of primitive cell',..
    hdr.pc_atom,'s',fid(1));

    PIL_print_mat('hdr.cc_vec, @f:f, lattice vectors of conventional cell',..
    hdr.cc_vec,'r',fid(1));
    PIL_print_mat('hdr.cc_sublat, @f:f, sublattice of conventional cell'..
    +' [index in pc(4),x, y, z, expan in cc(3)]',..
    hdr.cc_sublat,'r',fid(1),'off');
    PIL_print_mat('hdr.cc_atom, @f:f, atom type of conventional cell',..
    hdr.cc_atom,'s',fid(1));

    PIL_print_mat('hdr.slab_vec, @f:f, lattice vectors of slab cell',..
    hdr.slab_vec,'r',fid(1));
    PIL_print_mat('hdr.slab_sublat, @f:f, sublat of slab'..
    +' [sublat index in pc(1), x, y, z, proj on unit a3]',..
    hdr.slab_sublat,'r',fid(1),'off');  
    PIL_print_mat('hdr.slab_atom_type, @f:f, atom type of slab sublattice'..
    ,hdr.slab_atom,'s',fid(1))  

    PIL_print_mat('hdr.cc_unit, @f:f, pc in cc,'..
    +'[(b=1,n1,n2,n3),x,y,z,expan in cc]',hdr.cc_unit,'r',fid(1));

    PIL_print_mat('hdr.slab_unit, @f:f, pc in slab,'..
    +'[b=1,n1,n2,n3,x,y,z,expan in slab]',hdr.slab_unit,'r',fid(1));

    PIL_print_mat('hdr.coup_unit, @f:f, pc coupled to the R=0 slab,'..
    +'[(n1,n2,n3 in pc),(n1, n2, n3 in slab), uc_deg]',hdr.coup_unit,'i',fid(1));

    PIL_print_mat('hdr.state_info, @f:f, state info of slab WFs,'..
    +'[(wf label in pc), (n1, n2, n3 in pc)]',hdr.state_info,'i',fid(1));
    mclose(fid(1))
   
    if hdr.Verbosity=='more' then
        Hr_size=size(hdr.Hr_mat)
        printf('\n   => writing real part of Hr to '+..
        project_name+'_hdr_Hr_real.dat\n\n');
        fprintfMat(project_name+'_hdr_Hr_real.dat',..
        (real(full(hdr.Hr_mat)))'..
        ,'%f','real part of Hr, size=['..
        +string(Hr_size(1))+','+string(Hr_size(2))+']');
        
        printf('\n   => writing imag part of Hr to '+..
        project_name+'_hdr_Hr_imag.dat\n\n');
        fprintfMat(project_name+'_hdr_Hr_imag.dat',..
        (imag(full(hdr.Hr_mat)))'..
        ,'%f','imag part of Hr, size=['..
        +string(Hr_size(1))+','+string(Hr_size(2))+']');
    end
    
    // finishing program ===============================================
    save(project_name+'_hdr.sod','hdr');
    disp('{hdr}: finishing calculation ');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction
