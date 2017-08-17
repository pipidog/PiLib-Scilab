// **** Purpose ****
// PiLab Z2 invariant calculator
// **** variables ****
// << PiLab inputs >>
// [zti.Mesh]: 1x2, int
// <= mesh of the half-2D-BZ plane
// [zti.OccBand]: 1x1, int
// <= number of occupied bands
// [zti.Draw]: 1x1, string, 'on' / 'off'
// <= Whether to draw n-field, ''on'' / ''off''
// << PiLab outputs >>
// [zti.Z2_val]: 1x1 / 1x4, int
// => Z2 value, 3D: 1x4, 2D: 1x1
// [zti.Z2_plane]: 3x2 / NAN, int
// => planar Z2, 3D: [x0,x1; y0,y1;z0,z1], 2D: NAN
// [zti.n_field]: (total k of BZ) x (6) x (3 x 2)
// => n field configuration, [j1,j2,F, D2_A1, D1_A2, n_field]
//    n_field=F+D2_A1-D1_A2, if 2D, there is no (3x2) in size
// **** Version ****
// 11/24/2014 1st version
// 04/30/2015 add TR check before Z2 calculation
// 05/12/2015 change reload process, separate output parts
// **** Comment ****
// These formulas can be found in JPSJ 76 053702. Note that all of their
// formulas are based on full-periodic Bloch functions (FPBFs) rather
// than cell-periodic Bloch functions (CPBFs). In Comp Phys Comm 183
// 1849, they uses CPBFs, so periodic gauge has to be properly handled.
// Key component of this function is PIL_Z2_cal, check it for details.

function PiLab_zti(project_name)
    tic();
    disp('Starting {zti} calculation ...');
    disp('=========== Message ===========');

    // loading variables ===============================================
    disp('## loading variables ...');
    PiLab_loader(project_name,'zti','user','trim');
    load(project_name+'_zti.sod');
    load(project_name+'_lat.sod');
    load(project_name+'_hop.sod');
    load(project_name+'_scc.sod');

    // check variables =================================================
    disp('## checking variables ...')
    check_var=(length(zti.Mesh)==2)
    if check_var~=%t then
        disp('Error: PiLab_zti, zti.Mesh must be 1x2 integer vector!');
        abort;
    end
    check_var=(zti.OccBand <= length(hop.state_info(:,1)))
    if check_var~=%t then
        disp('Error: PiLab_zti, zti.OccBand must be less than total states');
        abort;
    end
    check_var=(zti.Draw=='on' | zti.Draw=='off')
    if check_var~=%t then
        disp('Error: PiLab_zti, zti.Draw must be ''on'' or ''off''');
        abort;
    end


    // Core part ========================================================
    disp('## running core part ...');
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
        zti.n_field=zeros(tot_k_plane,3,3,2);
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
        round([sum(zti.Z2_plane(1,:)),zti.Z2_plane(:,2)']),2); //[x0+x1,x1,x2,x3] mod 2
        disp('   -----------------');
        disp('   Total Z2_val= '+strcat(string(zti.Z2_val)));
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

    // output information ==============================================
    disp('## output information ...')
    if zti.Draw=='on'
        xsave(project_name+'_zti.scg');
        disp('   Output plot has saved to '+project_name+'_zti.scg'); 
    end

    select lat_size(1)
    case 2
        fid(1)=mopen(project_name+'_zti.plb','a+');
        PIL_print_mat('zti.Z2_val, @full, Z2 invariants'...
        ,zti.Z2_val,'i',fid(1)); 
        PIL_print_mat('zti.n_field, @full, n_field:'...
        +' [j1,j2,n_field]',zti.n_field,'r',fid(1));
        mclose(fid(1));
    case 3
        fid(1)=mopen(project_name+'_zti.plb','a+');
        PIL_print_mat('zti.Z2_val, @full, Z2 invariants'...
        ,zti.Z2_val,'i',fid(1)); 
        PIL_print_mat('zti.Z2_plane, @full, planar Z2'...
        +' [x0,x1;y0,y1;z0,z1]',zti.Z2_plane,'i',fid(1));
        for n=1:3
            for m=1:2
                select m
                case 1
                    print_fix='0' 
                case 2
                    print_fix='%pi' 
                end
                PIL_print_mat('zti.n_field(:,:,'+string(n)...
                +','+string(m)+'), @full, fields of G'+string(n)...
                +'@'+print_fix+': [j1,j2,n_field]'...
                ,zti.n_field(:,:,n,m),'r',fid(1));
            end
        end
        mclose(fid(1));
    end

    // finishing program ===============================================
    save(project_name+'_zti.sod','zti');
    disp('=============================');
    disp('Finishing {zti} calculation ...');
    disp('# time elapse= '+string(toc())+ ' seconds');
endfunction
