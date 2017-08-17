// **** Purpose ****
//  PiLab band state density analyzer(level 4)
// **** variables ****
// << PiLab inputs >>
// [ban_den.State]: Nx2, int, max: 4 states
// <= [k-point,state], the k-point label and eigenstate label in the ban 
//    output file. Can run N states simtansously.   
// [ban_den.Draw]: 1x1, string, 'on' / 'off'
// <= whether to draw results in bar  
// [ban_den.DrawType]: 1x1, string, 
// <= data to draw, 'site'/'orb'/'suborb'/'site_up'/'site_dn'
//    'orb_up','orb_dn'      
// << PiLab outputs>>
// [ban_den.site_den]: tot_site x N, real
// => density of each eigenstate on each site
// [ban_den.orb_den]: tot_orb x N, real
// => density of each eigenstate on each orbital
// [ban_den.suborb_den]: tot_suborb x N, real
// => density of each eigenstate on each basis states
// [ban_den.site_up_den]: tot_site x N, real
// => density of each eigenstate on each site w/ spin up
// [ban_den.site_dn_den]: tot_site x N, real
// => density of each eigenstate on each site w/ spin dn
// [ban_den.orb_up_den]: tot_orb x N, real
// => density of each eigenstate on each orbital w/ spin up
// [ban_den.orb_dn_den]: tot_orb x N, real
// => density of each eigenstate on each orbital w/ spin dn
// **** Version ****
// 05/10/2015 first built
// **** Comment ****

function PiLab_ban_den(project_name)
    tic();
    disp('Starting {ban_den} calculation ..');
    disp('=========== Message ===========');

    // loading variables ===============================================
    disp('## loading variables ..');
    PiLab_loader(project_name,'ban_den','user','trim');
    load(project_name+'_ban_den.sod'); 
    load(project_name+'_hop.sod');
    load(project_name+'_ban.sod');

    // check input =====================================================
    disp('## checking variables ..')
    check_var=length(ban_den.State(1,:))==2;
    if check_var~=%t then
        disp('Error: PiLab_ban_den, ban_den.State must be Nx2 !');
        abort;
    end
    check_var=length(ban_den.State(:,1))<=4;
    if check_var~=%t then
        disp('Error: PiLab_ban_den, ban_den.State must be less than 4 band states!');
        abort;
    end    
    check_var=(ban_den.Draw=='on' | ban_den.Draw=='off');
    if check_var~=%t then
        disp('Error: PiLab_ban_den, ban_den.Draw must be ''on'' or ''off'' !');
        abort;
    end
    check_var=(..
    ban_den.DrawType=='suborb' ..
    | ban_den.DrawType=='orb' ..
    | ban_den.DrawType=='site'..
    | ban_den.DrawType=='orb_up'..
    | ban_den.DrawType=='orb_dn'..
    | ban_den.DrawType=='site_up'..
    | ban_den.DrawType=='site_dn')
    if check_var~=%t then
        disp('Error: PiLab_ban_den,DrawType is not recognized !');
        abort;
    end    
    if length(grep(ban_den.DrawType,'_'))==1 then
        if (hop.Basis~='c') & (hop.Basis~='s') then
            disp('Error: PiLab_ban_den, hop.Basis must be ''c'' or ''s'''..
            +'to have spin component !');
            abort;
        end
    end    

    // Core part =======================================================
    disp('## running core part ..');

    state_info=hop.state_info; 
    state_info_text=hop.state_info_text; 
    clear hop

    tot_k=length(ban.k_band(1,:));
    tot_check=length(ban_den.State(:,1));

    tot_suborb=length(state_info(:,1));
    tot_orb=max(state_info(:,3));
    tot_site=max(state_info(:,2));

    ban_den.site_den=zeros(tot_site,tot_check);
    ban_den.orb_den=zeros(tot_orb,tot_check);
    ban_den.suborb_den=zeros(tot_suborb,tot_check);
    ban_den.site_up_den=zeros(tot_site,tot_check);
    ban_den.site_dn_den=zeros(tot_site,tot_check);
    ban_den.orb_up_den=zeros(tot_orb,tot_check);
    ban_den.orb_dn_den=zeros(tot_orb,tot_check);


    for n=1:tot_check
        disp('   running '+string(n)+'-th band states')
        ban_den.suborb_den(:,n)=(abs(ban.k_vec(:,ban_den.State(n,2)..
        ,ban_den.State(n,1)))).^2;
        for m=1:tot_suborb
            site=state_info(m,2);
            orb=state_info(m,3);
            if grep(state_info_text(m,5),',u')==1 then
                ban_den.site_up_den(site,n)=ban_den.site_up_den(site,n)..
                +ban_den.suborb_den(m,n);
                ban_den.orb_up_den(orb,n)=ban_den.orb_up_den(orb,n)..
                +ban_den.suborb_den(m,n);            
            elseif grep(state_info_text(m,5),',d')==1
                ban_den.site_dn_den(site,n)=ban_den.site_dn_den(site,n)..
                +ban_den.suborb_den(m,n);
                ban_den.orb_dn_den(orb,n)=ban_den.orb_dn_den(orb,n)..
                +ban_den.suborb_den(m,n); 
            end
            ban_den.site_den(site,n)=ban_den.site_den(site,n)..
            +ban_den.suborb_den(m,n);
            ban_den.orb_den(orb,n)=ban_den.orb_den(orb,n)..
            +ban_den.suborb_den(m,n);
        end
        //check results
        if abs(sum(ban_den.suborb_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.suborb_den='..
            +string(sum(ban_den.suborb_den(:,n))));
        end
        if abs(sum(ban_den.orb_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.orb_den='..
            +string(sum(ban_den.orb_den(:,n))));
        end
        if abs(sum(ban_den.site_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.site_den='..
            +string(sum(ban_den.site_den(:,n))));
        end
        if abs(sum(ban_den.site_up_den(:,n)+ban_den.site_dn_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.site_spin='..
            +string(sum(ban_den.site_up_den(:,n)+ban_den.site_dn_den(:,n))))
        end
        if abs(sum(ban_den.orb_up_den(:,n)+ban_den.orb_dn_den(:,n))-1) > 1e-5
            disp('Error: PiLab_ban_den, sum of ban_den.orb_spin='..
            +string(sum(ban_den.orb_up_den(:,n)+ban_den.orb_dn_den(:,n))));
        end
    end
    // Draw results  
    if ban_den.Draw=='on'
        for n=1:tot_check
            subplot(tot_check,1,n);
            select ban_den.DrawType
            case 'site'
                xlabel('site','fontsize',3);
                ylabel(string(n)+' den','fontsize',3);
                bar(ban_den.site_den(:,n));
            case 'orb'
                xlabel('orb','fontsize',3);
                ylabel(string(n)+' den','fontsize',3);
                bar(ban_den.orb_den(:,n));
            case 'suborb'
                xlabel('suborb','fontsize',3);
                ylabel(string(n)+' den','fontsize',3);
                bar(ban_den.suborb_den(:,n)); 
            case 'site_up'
                xlabel('site','fontsize',3);
                ylabel(string(n)+' up den','fontsize',3);
                bar(ban_den.site_up_den(:,n));
            case 'site_dn'
                xlabel('site','fontsize',3);
                ylabel(string(n)+' dn den','fontsize',3);
                bar(ban_den.site_dn_den(:,n));
            case 'orb_up'
                xlabel('orb','fontsize',3);
                ylabel(string(n)+' up den','fontsize',3);
                bar(ban_den.orb_up_den(:,n));
            case 'orb_dn'
                xlabel('orb','fontsize',3);
                ylabel(string(n)+' dn den','fontsize',3);
                bar(ban_den.orb_dn_den(:,n));
            end
            set(gcf(),'background',8)
        end

    end 

    // output ==========================================================
    disp('## output information ..')
    if ban_den.Draw=='on'
        xsave(project_name+'_ban_den_.scg');
        disp('   output plot has been saved to '+project_name+'_ban_den.scg');
    end

    // write to file
    fid(1)=mopen(project_name+'_ban_den.plb','a+');
    PIL_print_mat('ban_den.site_den, @full, the density on each site'..
    +'[den, eig_state]',ban_den.site_den,'r',fid(1));
    PIL_print_mat('ban_den.site_up_den, @full, the density on each site'..
    +' w/ spin up [den, eig_state]',ban_den.site_up_den,'r',fid(1));
    PIL_print_mat('ban_den.site_dn_den, @full, the density on each site'..
    +' w/ spin dn [den, eig_state]',ban_den.site_dn_den,'r',fid(1));
    PIL_print_mat('ban_den.orb_den, @full, the density on each orbital'..
    +' [den, eig_state]',ban_den.orb_den,'r',fid(1));
    PIL_print_mat('ban_den.orb_up_den, @full, the density on each orbital'..
    +' w/ spin up [den, eig_state]',ban_den.orb_up_den,'r',fid(1));
    PIL_print_mat('ban_den.orb_dn_den, @full, the density on each orbital'..
    +' w/ spin dn [den, eig_state]',ban_den.orb_dn_den,'r',fid(1));
    PIL_print_mat('ban_den.suborb_den, @full, the density on each suborbital'..
    +' [den, eig_state]',ban_den.suborb_den,'r',fid(1));
    mclose(fid(1));

    // finishing program ===============================================
    save(project_name+'_ban_den.sod','ban_den');
    disp('=============================');
    disp('Finishing {ban_den} calculation ..');
    disp('# time elapse= '+string(toc())+ ' seconds');
endfunction 
