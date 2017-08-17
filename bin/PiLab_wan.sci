// **** Purpose ****
// PiLab Wannier functions (Wannier90) interface (wan). 
// **** variables ****
//==== << PiLab inputs >> ====
//[wan.OutFile]: n x 1, string 
//<= file name of the wannier90 output of the wannier function
//[wan.HrFile]: n x 1, string
//<= file name of the wannier90 hr.dat file.
//[wan.Spin]: 1xn, string, 'r'/'s'/[]
//<= spin, 'r': relativistic, 's': scalar
//   PiLab will try to label the spin if possible. 
//[wan.Filter]: 1x1, real / []
//<= filter Hr_mat(:,:,n) smaller that given value, [] if no filter 
//[wan.Verbosity]: 1x1, string, 'less' / 'more'
//<= if 'less', H_r related information will not br printed.
//
//==== << PiLab outputs >> ====
//[wan.lat_vec]: 3x3, real
//=> lattice vectors
//[wan.rec_vec]: 3x3, real
//=> reciprocal lattice vectors
//[wan.atom_type]: nx1, string
//=> type of each atom
//[wan.sub_lat]: nx3, real
//=> sublattce positions
//[wan.wf_spatial]: nx4, ral
//=> [wannier function center (x,y,z), wannier spread]
//[wan.wf_site]: nx5,real
//=> nearest lattice site of each WFs. [b,n1,n2,n3,norm_err]
//[wan.state_info]: nx5, real
//=> [(x,y,z) in uc, spin]
//[wan.uc_index]: nx4, integer
//=> [uc_index(n1,n2,n3), degeneracy]
//[wan.Hr_mat]: tot_wf x tot_wf*tot_uc_index, sparse-real
//=> matrix of <n(0)|H|m(R)>, 
//   store in this way: Hr_mat(:,(n-1)*tot_wf+1:n*tot_wf)

//
//**** Version ****
// 01/19/2016 first built
// **** Comment ****
// This function combines level 1 ~ level 3 of PiLab (lat, hop, scc). 
// It directly reads the calculated resutls from W90. To run this, one
// must prepare the "wout" file and "hr.dat" file generated from 
// wannier90. This code will read them and generate them to PiLab 
// readable file formats. 

function PiLab_wan(project_name)
    disp('{wan}: starting calculation ');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{wan}: loading variables ');
    PiLab_loader(project_name,'wan','user','trim');
    load(project_name+'_wan.sod');

    // check variables  ================================================
    disp('{wan}: checking variables ')  
    if wan.Filter==[] | abs(wan.Filter) <1e-6 then
        disp('warning: wan.Filter reset to 1e-6');
        wan.Filter=1e-6;
    end
    check_var=(wan.Verbosity=='more' | wan.Verbosity=='less')
    if check_var~=%t then
        disp('Error: PiLab_wan, wan.Verbosity can be ''more'''..
        +' or ''less'' only!')
        abort
    end  
    disp('  => all variables paaed!')

    // core part =======================================================
    disp('{wan}: running core part ');
    // read wout file ------------------------------
    disp('  => reading Wannier functions')
    fid=mopen(wan.OutFile,'r');
    outfile_text=mgetl(fid,-1);

    // primitive vectors
    ind1=grep(outfile_text,'Lattice Vectors')
    wan.lat_vec=zeros(3,3)
    for n=1:3
        [t_len,t1,Vx,Vy,Vz]=msscanf(outfile_text(ind1(1)+n),'%s %f %f %f');
        wan.lat_vec(n,:)=[Vx,Vy,Vz];
    end

    //reciprical lattice
    ind1=grep(outfile_text,'Reciprocal-Space Vectors')
    wan.rec_vec=zeros(3,3)
    for n=1:3
        [t_len,t1,Vx,Vy,Vz]=msscanf(outfile_text(ind1(1)+n),'%s %f %f %f');
        wan.rec_vec(n,:)=[Vx,Vy,Vz];
    end

    //site coordinates
    ind1=grep(outfile_text,'Site       Fractional Coordinate')
    ind2=grep(outfile_text,'K-POINT GRID')
    tot_site=ind2(1)-ind1(1)-5
    wan.sub_lat=zeros(tot_site,3);
    wan.atom_type=[]
    for n=1:tot_site
        [t_len,t1,t2,d1,Vx1,Vy1,Vz1,t2,Vx2,Vy2,Vz2,t3]..
        =msscanf(outfile_text(ind1(1)+1+n),'%s %s %d %f %f %f %s %f %f %f %s');
        wan.atom_type(n)=t2;
        wan.sub_lat(n,:)=[Vx2,Vy2,Vz2];      
    end
    tot_sublat=length(wan.sub_lat(:,1));

    // read wannier function
    ind1=grep(outfile_text,'Final State')
    ind2=grep(outfile_text,'Final Spread')
    tot_wf=ind2(1)-ind1(1)-6
    wan.wf_spatial=zeros(tot_wf,4);
    for n=1:tot_wf
        [t_len,t1,t2,t3,t4,t5,t6,Vx1,t7,Vy1,t8,Vz1,t9,V2]..
        =msscanf(outfile_text(ind1(1)+n),'%2s %6s %3s %6s %f %1s %f %1s %f %1s %f %1s %f');
        wan.wf_spatial(n,:)=[Vx1,Vy1,Vz1,V2];      
    end
    mclose(fid)

    // find wf index and primitive sublattice correspondnece ----------
    disp('  => locating WFs nearest lattice sites')
    wan.wf_site=zeros(tot_wf,5);//[b,n1,n2,n3,norm_err]
    sublat_list_tmp=PIL_sublat_list(wan.lat_vec,wan.sub_lat,3);
    for n=1:tot_wf
        [wan.wf_site(n,1:4),wan.wf_site(n,5)]..
        =PIL_lat_index(wan.wf_spatial(n,1:3),wan.lat_vec,..
        wan.sub_lat,[],sublat_list_tmp);    
    end

    // read hr.dat -----------------------------------------------------
    disp('  => reading Hamiltonian')
    fid=mopen(wan.HrFile,'r');
    mgetl(fid,1);
    tot_wf=mfscanf(1,fid,'%f'); mgetl(fid,1);
    tot_uc=mfscanf(1,fid,'%f'); mgetl(fid,1);
    uc_deg=mfscanf(tot_uc,fid,'%f');mgetl(fid,1);

    hr_list=zeros(tot_uc,7);
    hr_list=mfscanf(tot_uc*tot_wf^2,fid,'%f %f %f %f %f %f %f\n')
    mclose(fid)
    // construct H(R) from Hmn(R)
    wan.uc_index=zeros(tot_uc,4);
    wan.Hr_mat=spzeros(tot_wf,tot_wf*tot_uc);

    uc_index_tmp=zeros(tot_uc,3);
    Hr_mat_tmp=zeros(tot_wf,tot_wf);
    index_2=0;
    count=0;
    printf('      total unit cell=%d\n',tot_uc)
    for n=1:tot_uc
        if pmodulo(n,100)==0
            printf('      %d/%d =%6.2f%% read\n',n,tot_uc,100*n/tot_uc)
        end
        // read uc_index and Hr
        index_1=index_2+1;
        index_2=index_1+tot_wf^2-1

        if hr_list(index_1,4)==1 & hr_list(index_1,5)==1 then
            uc_index_tmp=hr_list(index_1,1:3)
        else
            disp('Error!'); abort
        end  
        Hr_mat_tmp=clean(..
        matrix(hr_list(index_1:index_2,6)..
        +%i*hr_list(index_1:index_2,7),tot_wf,tot_wf)..
        ,wan.Filter);

        // eliminate useless all zero Hr
        if max(abs(Hr_mat_tmp))>1e-6 then
            count=count+1;
            wan.uc_index(count,:)=[uc_index_tmp,uc_deg(n)];
            wan.Hr_mat(:,(count-1)*tot_wf+1:count*tot_wf)=sparse(Hr_mat_tmp);
        end
    end
    wan.uc_index=wan.uc_index(1:count,:);
    wan.Hr_mat=wan.Hr_mat(:,1:count*tot_wf);

    // check spin state ------------------------------------------------
    disp('{wan}: checking spin states');
    wan.state_info=zeros(tot_wf,4);
    wan.state_info(:,1:3)=wan.wf_spatial(:,1:3);
    tot_state=length(wan.state_info(:,1));
    select wan.Spin
    case 'r'
        // generate TR_pair
        tmp=cat(2,[1:tot_state]',wan.wf_site(:,1:5));
        tmp=PIL_lsort(tmp,'c',[5,6,1:4],'i')
        TR_pair=gsort(matrix(tmp(:,1),2,-1)','r','i');
        // reorder state 
        tmp_index=find(TR_pair(:,1)>TR_pair(:,2));
        TR_pair(tmp_index,:)=flipdim(TR_pair(tmp_index,:),2);
        // check TR
        [TR_check,TR_mat]=PIL_TR_check(wan.lat_vec,wan.Hr_mat,..
        wan.uc_index(:,1:3),wan.uc_index(:,4),TR_pair,[1,2,3],1e-1);

        select TR_check
        case 1
            wan.state_info(TR_pair(:,1),4)=TR_pair(:,2);
            wan.state_info(TR_pair(:,2),4)=TR_pair(:,1);
            disp('  => TRS passed !');
        case 0
            disp('  => TRS failed !');
        end
    case 's'
        wan.state_info(1:tot_state/2,4)=-1;
        wan.state_info(1:tot_state/2,4)=+1;
    end

    // output information ==============================================
    disp('{wan}: output information ')
    // print suborblist
    fid=mopen(project_name+'_wan.plb','a+'); 
    PIL_print_mat('wan.lat_vec, @f:f, lattice vectors',..
    wan.lat_vec,'r',fid(1));
    PIL_print_mat('wan.rec_vec, @f:f, reciprocal lattice vectors',..
    wan.rec_vec,'r',fid(1));
    PIL_print_mat('wan.atom_type, @f:f, sublattice type',..
    wan.atom_type,'s',fid(1));
    PIL_print_mat('wan.sub_lat, @f:f, sublattice cartisan positions',..
    wan.sub_lat,'r',fid(1));
    PIL_print_mat('wan.wf_spatial, @f:f, [WF_center(x,y,z), WF_spread]',..
    wan.wf_spatial,'r',fid(1));
    PIL_print_mat('wan.wf_site, @f:f, [b,n1,n2,n3,norm_err]',..
    wan.wf_site,'r',fid(1),'off');
    PIL_print_mat('wan.state_info, @f:f, [(x,y,z) in uc, TR_pair]',..
    wan.state_info,'r',fid(1),'off');
    if wan.Verbosity=='more' then
        PIL_print_mat('wan.uc_index,  @f:f, [n1, n2, n3, degeneracy]',..
        wan.uc_index,'i',fid(1));
        for n=1:length(wan.uc_index(:,1))
            PIL_print_mat('wan.Hr_mat(:,(n-1)*tot_wf+1:n*tot_wf)'..
            +', @as:as, Hr matrix of '+string(n)+'-th uc_index',..
            wan.Hr_mat(:,(n-1)*tot_wf+1:n*tot_wf),'sp',fid);
        end
    end
    mclose(fid)

    // finishing program ===============================================
    save(project_name+'_wan.sod','wan');
    disp('{wan}: finishing calculation ');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction





