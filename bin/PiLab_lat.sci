// **** Purpose ****
// PiLab lattice generator 
// **** variables ****
// ==== << PiLab inputs >> ====
// [lat.LatVec]: 3x3/2x2/1x1, real
// <= primitive row vectors 
// [lat.Sublat]: nx3/nx2/nx1 , real
// <= sublattice position, 
// [lat.AtomType]: nx1, string or integer
// <= the atomic number of atomic name of each sublattice
//    e.g [73,73,33,33] or ['Ta','Ta','As','As']
// [lat.VecOrder]: 1x1, int
// <= order of primitive vector expansion for searching neighbors
// [lat.Tolerance]: 1x1, real
// <= criterion of identifying the order of neighbor. For bond length  
//    differnece small than this value will be considered as the same 
//    neighbor order. 
// [lat.Order]: 1x1, int
// <= select how many nn-orders will be print out.  
// ==== << PiLab outputs >> ====
// [lat.rec_vec]: 3x3/2x2/1x1, real
// => reciprocal lattice row vectors
// [lat.surr_site]: n x 9, real 
// => surrouding sites, [order, dist, sublatt, n1, n2, n3, x, y, z]
// **** Version ****
// 05/01/2014 1st version
// 05/12/2015 change reload process
// 02/09/2016 many variable names changed 
// 02/12/2016 lat.AtomType, lat.StrType added. support xsf output now.
// 04/06/2016 remove lat.Const to avoid confusion
// **** Comment ****
// 1. This code generates: 
//  1).reciprical lattice vectors (lat.rec_vec) 
//  2).information of surrounding sites (lat.surr_site)

function PiLab_lat(project_name)
    disp('{lat}: starting calculation ...');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{lat}: loading variables ...');
    // read current inputs
    PiLab_loader(project_name,'lat','user','trim');
    load(project_name+'_lat.sod');

    // check variables  ================================================ 
    disp('{lat}: checking variables ...') 
    // lat file
    check_var=(lat.VecOrder > 2);
    if check_var==%t then
        disp('Warning: PiLab_lat, lat.VecOrder is unusually large!');
    end
    check_var=(length(lat.LatVec(1,:))==3 & length(lat.Sublat(1,:))==3);
    if check_var~=%t then
        disp('Error: PiLab_lat, lat.LatVec and lat.Sublat must be 3D vectors');
        abort;
    end

    // core part ========================================================
    disp('{lat}: running core part ...');
    disp('  => generate structure file')

    atom_spec_tmp=[];
    for n=1:length(length(lat.AtomType))
        tmp=strsplit(lat.AtomType(n));
        tmp_idx=find(tmp=='*');
        tmp_rep=evstr(strcat(tmp(1:tmp_idx-1)));
        tmp_at=strcat(tmp(tmp_idx+1:$));
        atom_spec_tmp=cat(2,atom_spec_tmp,repmat(tmp_at,1,tmp_rep));
    end
    atom_spec=atom_spec_tmp;
    
    if length(length(atom_spec))~=length(lat.Sublat(:,1)) then
        disp('Error: PiLab_lat, number of lat.AtomType must be equal to'..
        +' total sublattice');
        abort;
    end
    
    PIL_crystal_xsf(project_name+'_lat',lat.LatVec,atom_spec..
    ,lat.Sublat)


    // generate reciprocal lattice vector
    disp('  => constructure reciprocal lattice')
    [lat.rec_vec]=PIL_recip_vec(lat.LatVec);
    // generate lat.surr_site
    disp('  => searching for neighbor sites')
    [lat.surr_site]=PIL_uc_nb(lat.LatVec...
    ,lat.Sublat,lat.VecOrder,lat.Tolerance,lat.Order);

    // output information ===============================================
    disp('{lat}: output information ...')
    fid(1)=mopen(project_name+'_lat.plb','a+'); 
    PIL_print_mat('lat.rec_vec, @f:f, the reciprocal lattice vectors'...
    ,lat.rec_vec,'r',fid(1));
    for n=1:length(lat.Sublat(:,1))
        PIL_print_mat('lat.surr_site('+string(n)+...
        '), @f:f, surrouding sites [order, dist, sublatt, n1, n2, n3, x, y, z]'...
        ,lat.surr_site(n),'r',fid(1),'off');
    end    
    mclose(fid(1));

    // finishing program ===============================================
    save(project_name+'_lat.sod','lat');
    disp('{lat}: finishing calculation ...');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction


