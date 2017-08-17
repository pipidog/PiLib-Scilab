// **** Purpose ****
// PiLab k-point information calculator. 
// This function is to calculate k-points using distributed parallel
// calculations, so dense k-mesh becomes possible. 
// **** Variables ****
//==== << PiLab inputs >> ====
// # System setup
//[kif.Command]: 1x1, string, 
//<= specify how do you submit a scilab job to the background in your 
//machine. The input file name and output file name should be replaced
//by @input_name and @output_name. PiLab will use correct name for them. 
//Note: In linux, it is usually:
//'scilab -nwni -f @input_name > @output_name &'
//In windows, it is usually:
//'start /b scilex -nwni -f @input_name > @output_name'

// # k-generation
//[kif.KptGen]: 1x1, string, 'grid'/'path'/'manual'
//<= Define how to generate k-points. If 'grid',PiLab will use 
//kif.SpanVec and kif.Mesh to generate grid k-points. 
//If 'path', PiLab use kif.Path and kif.Div to generate k-path
//If 'manual', PiLab will ignore kif.SpanVec and kif.Mesh and read the 
//k-points from a file called "kif_kpoint.dat". In this file, there 
//should be a nx3 matrix which lists all the k-points in PiLib variable 
//format.
//[kif.Format]: 1x1, string, 'red'/'cart'
//<= Define whether the k-points are in reduced or cartisian coordinate
// * If 'grid' --------------
//[kif.SpanVec]: 4x3, real
//<= Define 3 vectors to span a mesh box. The first row will be the
//origin and the next three are three positions from the origin
//to span the mesh box. Whether these points are defined in reduced
//or cartisna coordinate show be defined in kif.Format.
//[kif.Mesh]: 1x3, int
//<= Define mesh number of each vector. Note that, the mesh will 
//include head and tail of each vector. If any mesh number=1, means
//there is no mesh. If mesh number=2, means only the origin and the
//end point are included in that direction.     
// * If 'path' ---------------
//[kif.Path]: nx3, real
//<= k_path in reduced coordinates
//[kif.Div]: 1x1, int
//<= k_mesh along k_path of each unit length k

// # Job allocation   
//[kif.Allocation]: 1x1, int
//<= how many allocations to disturbut this calculation.
//[kif.Thread]: 1x1, int
//<= how many threads to run at the same time
//[kif.MaxRAM]: 1x1, real
//<= set the Max RAM usage for the calculation. If the PiLab estimeted
//RAM is larger than this, calculation will be stopped in the beginning. 
//If so, try to increase Kif.Allocation so each allocation will need
//less RAM. 
//==== << PiLab output >> ====
//[kif.k_point]: nx3, real
//=> k-points
//[kif.k_num]: 1 x tot_alloc, int
//=> number of k-points in each allocation
// **** Version ****
// Apr 29, 2016
// **** Comment ****

function PiLab_kif(project_name)
    disp('{kif}: starting calculation ');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{kif}: loading variables ');
    PiLab_loader(project_name,'kif','user','trim');
    load(project_name+'_kif.sod');
    load(project_name+'_ham.sod');
    disp('  => all variables loaded');

    // variable check ===================================================
    disp('{kif}: checking variables ')
    if grep(kif.Command,'&')==[] & grep(kif.Command,'start /b')==[] then
        disp('Warning: PiLab_kif, kif.Command must have send-to-background'..
        +' command, e.g. & (linux) / start /b (windows) !');
    end
    if kif.KptGen~='grid' & kif.KptGen~='path' & kif.KptGen~='manual' then
        disp('Error: PiLab_kif, kif.KptGen can only be ''grid/path/manual''');
        abort
    end
    if kif.Format~='red' & kif.Format~='cart' then
        disp('Error: PiLab_kif, kif.Format can only be ''red'' or ''cart''');
        abort
    end
    select kif.KptGen
    case 'grid'
        if length(kif.SpanVec(1,:))~=3 | length(kif.SpanVec(:,1))~=4 then
            disp('Error: PiLab_kif,kif.SpanVec must be 4x3 matrix');
            abort
        end
        if length(kif.Mesh)~=3 | find(kif.Mesh<=0)~=[] then
            disp('Error: PiLab_kif,kif.Mesh be a 1x3 and all >=1');
            abort
        end
    case 'path'
        if length(kif.Path(1,:))~=3 then
            disp('Error: PiLab_kif, kif.Path must be nx3 !');
            abort
        end
        if length(kif.Div)~=1 | kif.Div-fix(kif.Div)~=0 then
            disp('Error: PiLab_kif, kif.Div must be an integer !');
            abort
        end
    end
    disp('  => all variables passed')

    // core Part ========================================================
    // generate k-point
    disp('{kif}: running core part ');
    disp('  => generating k-points');
    select kif.KptGen
    case 'grid'
        [kif.k_point]=PIL_k_mesh(ham.rec_vec,kif.Mesh,kif.SpanVec,kif.Format,'on');
    case 'path'
        if kif.Format=='red' then
            [kif.k_point,kif.k_div]=PIL_k_path(kif.Path*ham.rec_vec,kif.Div,'unit');
        else
            [kif.k_point,kif.k_div]=PIL_k_path(kif.Path,kif.Div,'unit');
        end
    case 'manual'
        fid=mopen('kif_kpoint.dat','r');
        kif.k_point=PIL_read_mat(fid);
        mclose(fid);
        if kif.Format=='red' then
            kif.k_point=kif.k_point*ham.rec_vec;
        end
    end

    // estimate memory usage for each allocation
    disp('  => estimate memory usage ');
    tot_k=length(kif.k_point(:,1));
    tot_state=length(ham.state_info(:,1));
    n_k=ceil(tot_k/kif.Allocation);
    tot_ele=(5*tot_state*tot_state+tot_state)*n_k; // Vk+Hk part
    // estimated RAM in MB=matrix in RAM + file size when store+10% buffer
    est_mem=((tot_ele/1e7)*76+(tot_ele/1e6)*15.2)*1.1; 
    disp('     estimated memory per allocation= '+string(est_mem)+'MB');
    if est_mem > 2040.0 then
        disp('Warning: estimated memory is larger than 2GB !');
        disp('         Reset to 2GB, job could fail!');
        est_mem=2040.0;
    end
    if est_mem*kif.Thread > kif.MaxRAM then
        disp('Error: PiLab_kif, estimeted memory is larger than kif.MaxRAM!');
        disp('       Please increase kif.Allocation.');
        abort;
    else
        disp('     total memory needed is '+string(est_mem*kif.Thread)+'MB'); 
    end
    
    disp('  => estimate file size ')
    est_size=(n_k*tot_state^2)/65000;
    disp('     estimated Vk size per allocation= '+string(est_size)+'MB')
    if  est_size > 350.0 then
        disp('Error: PiLab_kif, estimated file size > 350.0 MB !');
        disp('       Please increase allocations!');
        abort
    end
    clear ham

    // perform distrubited parallel calculation
    disp('  => running distributed calculation ');
    disp('     check ./kif/runtime.log for runtime status');
    [kif.k_alloc]=PIL_Hk_solver_dist(pwd(),project_name,'kif',kif.Command,..
    kif.k_point,3,kif.Allocation,kif.Thread,est_mem);

    // output information ==============================================
    disp('{kif}: output information ');
    fid(1)=mopen(project_name+'_kif.plb','a+');
    PIL_print_mat('kif.k_alloc, @f:f, number of k-point at each allocation'..
    +' [tot_k,ini_k,final_k]',kif.k_alloc,'i',fid(1));
    if kif.KptGen=='path' then
        PIL_print_mat('kif.k_div, @f:f, k_div',kif.k_div,'i',fid(1));
    end
    PIL_print_mat('kif.k_point, @f:f, k_points',kif.k_point,'r',fid(1));
    mclose(fid(1));

    // finishing program ===============================================

    save(project_name+'_kif.sod','kif');
    disp('{kif}: finishing {kif} calculation ');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction

