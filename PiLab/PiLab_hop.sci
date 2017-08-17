// **** Purpose ****
// PiLab hopping generator (level 2)
// **** variables ****
// ==== << PiLab inputs >> ====
// [hop.SiteOrb]: n x 2, integer 
// <= specify orbital of each site, nx2, [site, l]
// [hop.Order]: 1x1 integer
// <= order of nearest coupling 
// [hop.SKint]: n x 7, real 
// <= SK parameters, [Orb1,Orb2,nn_order,ts,tp,td,tf]
// [hop.LS]: 1xn, real
// <= strength of LS coupling on each SiteOrb. 
// [hop.Filter]: 1x1, real
// <= filter of small hopping elements, 
// [hop.Basis]: 1x1, string
// <= basis of Hamiltonian, 'c', 's', 'rc', 'rs'
// [hop.SelState]: 1xn, integer
// <= select states by inputting their state label 
// [hop.OnsiteE]: 1xn,, real 
// <= Onsite energy of selected states (given by their order)
// ==== << PiLab outputs >> ====
// [hop.state_info_text]: total state x 5, string  
// => state info in text format
// phop.state_info]: total state x 5, int 
// => state info in num format, 'i'
// [hop.LS_mat]: n x 3, t-sp real
// => LS coupling matrix
// [hop.onsite_E]: n x 3, t-sp real
// => onsite energy matrix
// [hop.hop_size]: total sublattice x 4, int
// => size of hop_mat
// [hop.hop_mat]: list(n) x (n,m,p), a-sp 
// => hopping matrix 
// **** Version ****
// 05/01/2014 1st version
// 05/12/2015 change reload process
// **** Comment ****
// 1. This code generates: 
//  1). state information in text format (hop.state_info_text)
//  2). state information in value format (hop.state_info)
//  3). LS coupling matrix (hop.LS_mat, sparse)
//  4). the size of hopping matrix (hop.mat_size)
//  5). the hopping matrix (hop.hop_mat, sparse)

function PiLab_hop(project_name)
    disp('{hop}: starting calculation ...');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{hop}: loading variables ...');
    PiLab_loader(project_name,'hop','user','trim');
    load(project_name+'_hop.sod');
    load(project_name+'_lat.sod');
    
    // check variables  ================================================
    disp('{hop}: checking variables ...')    
    check_var=(hop.Order <=lat.Order);
    if check_var~=%t then
        disp('Error: PiLab_hop, hop.Order must <= lat.Order!');
        abort;
    end
    check_var=(find(hop.SiteOrb(:,1)>length(lat.Sublat(:,1)))==[]);
    if check_var~=%t then
        disp('Error: PiLab_hop, hop.SiteOrb has wrong site index!');
        abort;
    end
    for n=1:length(lat.Sublat(:,1))
        check_var=(length(find(hop.SiteOrb(:,1)==n))~=0);
        if check_var~=%t then
            disp('Error: PiLab_hop, hop.SiteOrb has unassigned site!');
            abort;
        end
    end
    check_var=(length(hop.SKint(1,:))==7);
    if check_var~=%t then
        disp('Error: PiLab_hop, hop.SKint should has 7 '...
        +'index for a single coupling!');
        abort;
    end
    check_var=(max(hop.SKint(:,1:2))<=length(hop.SiteOrb(:,1)));
    if check_var~=%t then
        disp('Error: PiLab_hop, hop.SKint has wrong identifiers!');
        abort;
    end
    check_var=(length(hop.LS)==length(hop.SiteOrb(:,1)));
    if check_var~=%t then
        disp('Error: PiLab_hop, hop.LS must have the same '..
        +' length as hop.SiteOrb(:,1))');
        abort
    end
    check_var=(hop.Basis=='s' | hop.Basis=='c' ...
    | hop.Basis=='rs' | hop.Basis=='rc');
    if check_var~=%t then
        disp('Error: PiLab_hop, hop.Basis can only be '...
        +'''s'', ''c'', ''rs'', ''rc''!');
        abort;
    end
    check_var=(length(hop.SelState)==length(hop.OnsiteE));
    if check_var~=%t then
        disp('Error: PiLab_hop, hop.SelState and '...
        +'hop.OnsiteE have incosistent dimension!');
        abort;
    end
    
    // core part ========================================================
    disp('{hop}: running core part ...');
    disp('  => generating state information')
    // enforce orbitab identifier used in PIL_hop_mat
    hop.SiteOrb=hop.SiteOrb(:,[1,1,2]);  
    // to be the input order
    hop.SiteOrb(:,2)=[1:length(hop.SiteOrb(:,2))]';     

    // generate hopping integrals --------------------------------------
    disp('  => generating Slaster-Koster hopping integrals')
    hop.SiteOrb=gsort(hop.SiteOrb,'lr','i');
    [hop.state_info,hop.hop_mat]=PIL_hop_mat(lat.surr_site,hop.SiteOrb...
    ,hop.SKint,hop.Basis,hop.Order,hop.Filter);

    //generate LS coupling ---------------------------------
    disp('  => generating spin-orbit coupling')
    hop.LS_mat=[];
    for n=1:length(hop.SiteOrb(:,1))
        select hop.SiteOrb(n,3)
        case 0
            hop.LS_mat=PIL_dirsum(hop.LS_mat...
            ,PIL_LS_coup('s',hop.LS(hop.SiteOrb(n,2)),hop.Basis));
        case 1
            hop.LS_mat=PIL_dirsum(hop.LS_mat...
            ,PIL_LS_coup('p',hop.LS(hop.SiteOrb(n,2)),hop.Basis));
        case 2
            hop.LS_mat=PIL_dirsum(hop.LS_mat...
            ,PIL_LS_coup('d',hop.LS(hop.SiteOrb(n,2)),hop.Basis));
        case 3
            hop.LS_mat=PIL_dirsum(hop.LS_mat...
            ,PIL_LS_coup('f',hop.LS(hop.SiteOrb(n,2)),hop.Basis));
        end
    end

    // generate  onsite energy ---------------
    disp('  => generating onsite energy')
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

    // output information ==============================================
    disp('{hop}: output information ...')
    // print suborblist
    fid(1)=mopen(project_name+'_hop.plb','a+'); 
    PIL_print_mat('hop.state_info_text, @f:f'...
    +', [state_label, site, identifier, l, SubOrb_text] '...
    ,hop.state_info_text,'s',fid(1));

    PIL_print_mat('hop.state_info, @f:f, [state_label'...
    +', site, identifier, l, SubOrb] ',hop.state_info,'i',fid(1));

    PIL_print_mat('hop.LS_mat, @f:ts, LS coupling matrix',...
    sparse(hop.LS_mat),'sp',fid(1));

    PIL_print_mat('hop.onsite_E, @f:ts, onsite energy matrix'...
    ,sparse(hop.onsite_E),'sp',fid(1)); 

    PIL_print_mat('hop.hop_size, @f:f, size of hop.hop_mat'...
    +', [sublatt, size(hop.hop_mat(n))]',hop.hop_size,'i',fid(1));

    for n=1:size(lat.surr_site)
        for m=1:length(lat.surr_site(n)(:,1))-1
            PIL_print_mat('hop.hop_mat('+string(n)+')(:,:,'+string(m)...
            +'), @f:as, hop_mat between site-'...
            +string(n)+' and its '+string(m)+'-th neighbor'...
            ,sparse(squeeze(hop.hop_mat(n)(:,:,m))),'sp',fid(1)); 
        end
    end
    mclose(fid(1));

    // finishing program ===============================================
    save(project_name+'_hop.sod','hop');
    disp('{hop}: finishing calculation ...');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction


