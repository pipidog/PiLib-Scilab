// **** Purpose ****
// PiLab Floquet Operator generator
// **** variables ****
// ==== << PiLab inputs >> ====
// [flq.Frequency]: 1x1 real
// <= field frequency, 
// [flq.Order]: 1x1, int
// <= order of photon process, 
// [flq.Amplitude]: 1x1 / 1x2 / 1x3, real
// <= AC amplitude,  
// [flq.Phase]: 1x1 / 1x2 / 1x3, real
// <=AC phase
// ==== << PiLab outputs >> ====
// [flq.state_info]: tot flq state x 6, int
// => [state_label, order, site, identifier, l, SubOrb]
// [flq.H_local]: list(tot order) x n x 3, t-sp
// => onsite energy of a particualr order 
// [flq.hop_size]: (tot_order*tot_sublatt) x 5
// => size of flq.hop_mat, [order+1,sublatt,hop_mat_size]
// [flq.hop_mat]: list(tot order)(tot_sublatt) x hop_mat(:,:,:), a-sp
// => floquet hopping matrix 
// **** Version ****
// 07/03/2014 1st version
// 05/12/2015 change reload process
// 06/04/2015 add flq.E_Fermi , flq.E_gap
// **** Comment ****
// 1. In Floquet theory, the hopping integral will be renormalized 
//    by the Floquet integral. This code calculates the renormalized
//    hopping for assgined order of photon process. Once you get the
//    renormalized hopping integrals, you can construct and diagonalize
//    the Floquet operator on different k-point by using PiLab_flq_ban
//    or PiLab_flq_dsa, etc. 
// 2. The code requires the external periodic ptential is homogeneous. 
//    for 2D and 1D cases, it is easy to be achieved by an EM wave with
//    assumption of z=0. For 3D, EM wave doesn't work becaues it has a
//    propagation direction. 
// 3. The form of vector potential
//    A=[Ax*cos(w*t+s1),Ay*cos(wt+s2),Az*cos(wt+s3)], r=[rx,ry,rz]
//    s1,s2,s3: the phase terms
// 4. For 2D case (Az=0) with Ax=Ay, then
//    (s2-s1)=0 --> linear (right up to left down)
//    (s2-s1)=%pi/2 --> circular
//    (s2-s1)=%pi --> linear (left up to right down)  

function PiLab_flq(project_name)
    disp('{flq}: starting calculation ...');
    c1=clock();
    printf('\n');
    printf('   => start time: %4d/%02d/%02d %02d:%02d:%02d\n',c1);
    // loading variables ===============================================
    disp('{flq}: loading variables ...');
    PiLab_loader(project_name,'flq','user','trim');
    load(project_name+'_flq.sod');
    load(project_name+'_lat.sod');
    load(project_name+'_hop.sod');

    // check variable ==================================================
    disp('{flq}: checking variables ...')
    check_var=(length(flq.Frequency)==1);
    if check_var~=%t then
        disp('Error: PiLab_flq, flq.Frequency must be 1x1 real !');
        abort;
    end
    check_var=(length(flq.Order)==1 & (fix(flq.Order)-flq.Order)==0 );
    if check_var~=%t then
        disp('Error: PiLab_flq, flq.Order must be 1x1 integer !');
        abort;
    end
    check_var=(length(flq.Amplitude)==length(lat.LatVec(:,1)));
    if check_var~=%t then
        disp('Error: PiLab_flq, flq.Amplitude must '...
        +'have the same dimension with your system !');
        abort;
    end
    check_var=(length(flq.Phase)==length(lat.LatVec(:,1)));
    if check_var~=%t then
        disp('Error: PiLab_flq, flq.Phase must have the '...
        +'same dimension with your system !');
        abort;
    end
    
    // Core Part =======================================================
    disp('{flq}: running core part ...');
    // calculate renormalized matrix elements -------------------------
    // flq.H_local(order_index)
    // flq.hop_mat(order_index)(sublatt_index)
    disp('  => generating Floquet hopping integrals')
    flq.H_local=list()
    flq.hop_mat=list();
    for p=1:flq.Order+1 // run photon process order, 0~flq.Order
        // Oniste renormalization
        flq.H_local(p)=hop.onsite_E+hop.LS_mat;
        flq.H_local(p)=flq.H_local(p)...
        *PIL_flqint(flq.Frequency,p-1,flq.Amplitude,flq.Phase,[0,0,0]);

        // hopping renormalization
        flq.hop_mat(p)=list();
        for n=1:length(hop.hop_size(:,1)) // run all sublattice
            r1=lat.surr_site(n)(1,7:9);
            flq.hop_mat(p)(n)=hop.hop_mat(n);
            for m=1:hop.hop_size(n,4) // run all nearest neighbor
                // renormalize hopping integrals
                r2=lat.surr_site(n)(m+1,7:9);    
                flq.hop_mat(p)(n)(:,:,m)=flq.hop_mat(p)(n)(:,:,m)...
                *PIL_flqint(flq.Frequency,p-1...
                ,flq.Amplitude,flq.Phase,r2-r1);
            end
        end
    end

    // generate Floquet state index ------------------------------------
    disp('  => generating Floquet state information')
    flq.state_info=[];
    for n=-flq.Order:flq.Order
        flq.state_info=cat(1,flq.state_info...
        ,cat(2,n*ones(hop.state_info(:,1)),hop.state_info(:,2:$)));
    end
    flq.state_info=cat(2,[1:length(flq.state_info(:,1))]'...
    ,flq.state_info);
    // state info in text
    // generate flq.state_info_text
    suborb_list=PIL_suborb_list(hop.Basis);
    flq.state_info_text=cat(2,string(flq.state_info(:,1:5)),...
    suborb_list(flq.state_info(:,$)))

    // generate size of flq.hop_mat for PiLab_loader -------------------
    flq.hop_size=[];
    for p=1:flq.Order+1
        flq.hop_size=cat(1,flq.hop_size,cat(2,...
        p*ones(length(hop.hop_size(:,1)),1),hop.hop_size));
    end

    // output information ==============================================
    disp('{flq}: output information ...')
    fid=mopen(project_name+'_flq.plb','a+');
    // Floquet state text
    PIL_print_mat('flq.state_info_text, @f:f, [state_label'...
    +', order, site, identifier, l, SubOrb_text]'...
    ,flq.state_info_text,'s',fid);

    // Floquet state labels
    PIL_print_mat('flq.state_info, @f:f, [state_label'...
    +', order, site, identifier, l, SubOrb]'...
    ,flq.state_info,'i',fid);
    
    // print flq.H_local 
    for p=1:flq.Order+1
        PIL_print_mat('flq.H_local('+string(p)...
        +'), @f:ts, Floquet H_onsite of order '...
        +string(p-1), sparse(flq.H_local(p)),'sp',fid);
    end
    // print flq.hop_mat and its size
    PIL_print_mat('flq.hop_size, @f:f, size of flq.hop_mat'...
    +' [order+1,sublatt,hop_mat_size]'...
    ,flq.hop_size,'i',fid);
    for p=1:flq.Order+1
        for n=1:length(hop.hop_size(:,1))
            for m=1:hop.hop_size(n,4)
                PIL_print_mat('flq.hop_mat('+string(p)...
                +')('+string(n)+')(:,:,'...
                +string(m)+'), @f:as, Floquet hop_mat('...
                +string(n)+')(:,:,'+string(m)+') of order '+string(p-1)...
                ,sparse(squeeze(flq.hop_mat(p)(n)(:,:,m))),'sp',fid);
            end
        end
    end
    mclose(fid);

    // finishing program ===============================================
    save(project_name+'_flq.sod','flq');
    disp('{flq}: finishing calculation ...');
    disp('  => time elapse '+string(etime(clock(),c1))+ ' seconds');
endfunction
