// **** Purpose ****
// PiLab Chern number caculator
// **** variables ****
// << PiLab inputs >>
// [chn.Mesh]: 1x2, int
// <= k-space mesh
// [chn.OccBand]: 1x1, int
// <= number of occupied bands
// [chn.Kdiff]: 1x2, real
// <= Differential vector to avoid divergence
// << PiLab Outputs >>
// [chn.tot_Chern]: 1x1, int
// => total Chern number
// [chn.ban_Chern]: OccBand x 1, int
// => Chern number of each occupied band
// [chn.Fk_field]: occBand x total k of BZ
// => the Chern field of each band at each k-point
// **** Version ****
// 08/31/2014 1st version
// 12/02/2014 separate Chern calculation to independent functions
//            totally rewrite to highly improve effiency
// 05/12/2015 change reload process
// **** Comment ****
// see JPSJ 74.1674 (2005), eq.(6)~ eq.(9)

function PiLab_chn(project_name)
    tic();
    disp('Starting {chn} calculation ...');
    disp('=========== Message ===========');
    
    // loading variables ===============================================
    disp('## loading variables ...');
    PiLab_loader(project_name,'chn','user','trim');
    load(project_name+'_chn.sod');
    load(project_name+'_lat.sod');
    load(project_name+'_hop.sod');
    load(project_name+'_scc.sod');
    
    // check variables =================================================
    disp('## checking variables ...')
    check_var=(length(lat.Primitive(1,:))==2);
    if check_var~=%t then
        disp('Error: PiLab_chn, Chern numbers only survive in 2D');
        abort;
    end
    
    check_var=(length(chn.Mesh)==2 & find(chn.Mesh<=0)==[]);
    if check_var~=%t then
        disp('Error: PiLab_chn, chn.Mesh has wrong '...
        +'dimension or non-positive integers!');
        abort;
    end
    
    check_var=((chn.OccBand >0) & (chn.OccBand-fix(chn.OccBand))==0)
    if check_var~=%t then
        disp('Error: PiLab_chn, chn.OccBand must be a positive integer!');
        abort;
    end
    
    check_var=(length(chn.Kdiff)==2)
    if check_var~=%t then
        disp('Error: PiLab_chn, chn.Kdiff has wrong dimension!');
        abort;
    end
    
    // Core Part =======================================================
    disp('## running core part ...');
    [chn.k_point,chn.Chern_val,chn.Fk_field]...
    =PIL_Chern_cal(lat,hop,scc,[],chn.Mesh,chn.OccBand,chn.Kdiff)

    // output calculation results ======================================
    disp('## output information ...')
    fid(1)=mopen(project_name+'_chn.plb','a+');
    PIL_print_mat('chn.Chern_val, @full, total Chern number below Ef'...
    ,chn.Chern_val,'r',fid(1)); ;
    PIL_print_mat('chn.k_point, @full, k-point inside BZ'...
    ,chn.k_point,'r',fid(1));
    PIL_print_mat('chn.Fk_field, @full, F-field at each'...
    +' k-point [F(k)/(2*%pi*%i)]',chn.Fk_field,'r',fid(1));
    mclose(fid(1));
    disp('Total Chern number= '+string(clean(chn.Chern_val)));
    
    // finishing program ===============================================
    save(project_name+'_chn.sod','chn');
    disp('=============================');
    disp('Finishing {chn} calculation ...');
    disp('# time elapse= '+string(toc())+ ' seconds');
endfunction
