// **** Purpose ****
// plot fat band structure
// **** Variables ****
// [Ek]: real, tot_ban x tot_k
// <= eigenvalues of the band structure, must arragned as tot_ban x tot_k
// [Ek_weight]: real, tot_ban x tot_k x tot_proj
// <= projected weight on each eigenvalues 
// [k_div]: int, 1 x tot_div
// <= location of k where a divider should be placed, e.g: [1,35, 61,82,100]
//      first and last should be at the boundaries
// [E_bound]: int, 1x2
// <= E_bound for the plot, e.g. [-5,5]
// [state_grp]: list, optional, default:[]
// <= how to group states in plot. type of list. 
//    e.g: list([1:4],[7:13]), then this function will combine the weight
//         of state 1~4 as first plot and 7:13 as second plot. 
//    if [], don't group any state.
// [state_info]:strings, tot_state x 1, optional, default=[]
// <= information of each state. 
// [k_label]: length(k_div)x1, string
// <= labels of each high symmetry points,e.g. ['$\Gamma$','$\Delta$','L'] 
// [desc]: str, 1x1, optional, default=''
// <= desription of the plots. will be used as the title of each sub figures. 
// [ini_fig_num]: int, 1x1, optional, default=0
// <= the initial figure ID for the series of plot. it is to avoid replacement of your previous plots
//       if you call this function multiple times in a code. 
// [marker_size]: int, 1x1, optiona, default=6
// <= size of the marter
// [subplot_size]: int, 1x2, optional, default=[]
// <= size of the subplot. if [], won't use subplot
// [c_map]: real, 2x3, values: 0~1, optional, default=[]
// <= color map, if [], use default, i.e. blue and gray
// [w_range]: real, 1x1, optional, default=[] 
// <= renomarlization of the wieght
// **** Version ****
// Apr 17, 2016 first built
// 06/03/2016 add cmax parameter
// **** Comment ****

function PIL_fatband_plot(Ek,Ek_weight,Ef,k_div,E_bound,state_grp,state_info,k_label,desc,ini_fig_num,marker_size,subplot_size,c_map,w_range,k_norm)
    default_val=['state_grp=[]','state_info=[]','k_label=[]','desc=[]',..
    'ini_fig_num=[]','marker_size=[]','subplot_size=[]',..
    'c_map=[]','w_range=[]','k_norm=''off'''];
    [lhs,rhs]=argn();
    if (rhs>5) &(rhs<15) then
        for n=rhs+1:15
            execstr(default_val(n-5))
        end
    end
    tot_k=length(Ek(1,:));
    tot_ban=length(Ek(:,1));

    // add ignored input variable
    if k_label==[] then
        k_label='k'+string(linspace(1,length(k_div),length(k_div)));
    end
    if desc==[] then
        desc='Title'
    end
    if ini_fig_num==[] then
        ini_fig_num=0
    end
    if marker_size==[] then
        marker_size=6
    end
    if k_norm=='on' then
        k_range=linspace(0,1,tot_k)
    else
        k_range=linspace(1,tot_k,tot_k)
    end

    if state_info==[] then
        Ek_weight_size=size(Ek_weight);
        if length(Ek_weight_size)==3 then
            state_info=string([1:Ek_weight_size(3)]');
        else
            state_info='1'
        end
    end

    if length(length(k_label))~=length(k_div) then
        disp('Error: PIL_fatband_plot, size(k_label)~=size(k_div)!');
        abort 
    end
    // combine Ek_weight group
    if state_grp~=[] then
        Ek_weight_grp=zeros(length(Ek(:,1)),length(Ek(1,:)),..
        length(state_grp));
        printf('\n');
        for n=1:length(state_grp)
            printf('projected states on plot-%d:\n',n);
            printf('   %s\n',state_info(state_grp(n)))
            for m=1:length(state_grp(n))
                Ek_weight_grp(:,:,n)=..
                Ek_weight_grp(:,:,n)+Ek_weight(:,:,state_grp(n)(m))
            end
        end
        Ek_weight=Ek_weight_grp;
        clear Ek_weight_grp;
    else
        for n=1:length(state_grp)
            printf('projected states on plot-%d:\n',n);
            printf('   %s\n',string(n));
        end
    end



    // pick bands in the window
    if E_bound~=[] then
        sel_ban=find(Ek(:,1)>=E_bound(1) & Ek(:,1)<=E_bound(2));
        if sel_ban==[] then
            HOMO_ban=max(find(Ek(:,1)<0));
            sel_ban=HOMO_ban-4:HOMO_ban+5;
        else
            if sel_ban(1)>=6 then
                sel_ban=[sel_ban(1)-5:sel_ban(1)-1,sel_ban];
            else
                sel_ban=[1:sel_ban($)];
            end

            if tot_ban-sel_ban($)>=5 then
                sel_ban=sel_ban(1):sel_ban($)+5;
            else
                sel_ban=sel_ban(1):tot_ban;
            end
        end

    else
        sel_ban=1:tot_ban;
    end
    tot_sel_ban=length(sel_ban);

    // plot band structure
    tot_proj_task=size(Ek_weight(1,1,:));
    if length(tot_proj_task)==3 then
        tot_proj_task=tot_proj_task(3);
    else 
        tot_proj_task=1;
    end

    if c_map~=[];
        c_map=flipdim(PIL_k_path([c_map(2,:);c_map(1,:)],32),1);
    else
        c_map_low=[0.75 0.75 0.75]
        c_map_high=[0.0 0.0 0.6]
        c_map=flipdim(PIL_k_path([c_map_high;c_map_low],32),1);
    end     

    printf('\n')
    for n=0:tot_proj_task        
        if n==0 then  // original band
            scf(ini_fig_num+0);
            plot(k_range',Ek','b','thickness',3);
            font_size=4;
            linewidth=4
            set(gcf(),'background',8);

        else  // projected band
            if subplot_size~=[] then
                scf(ini_fig_num+1);
                subplot(subplot_size(1),subplot_size(2),n)
                font_size=4;
                linewidth=2;
            else
                scf(ini_fig_num+n);
                font_size=4;
                linewidth=2;
            end       
            // scatter plot
            if grep(getversion(),'scilab-6')~=[] then
                // new version for scilab 6.x
                xset("colormap",c_map);                
                x=matrix(repmat(k_range,tot_sel_ban,1),-1,1)
                y=matrix(Ek(sel_ban,:),-1,1)
                w=matrix(Ek_weight(sel_ban,:,n),-1,1);
                [w_new,w_sort]=gsort(w,'g','i')
                scatter(x(w_sort),y(w_sort),marker_size,w_new,'fill')
                colorbar(min(w),max(w));
            else
                // old version for scilab 5.x
                PIL_scatter_plot(..
                matrix(repmat(k_range,tot_sel_ban,1),-1,1),..
                matrix(Ek(sel_ban,:),-1,1),matrix(Ek_weight(sel_ban,:,n),-1,1),..
                marker_size,c_map,w_range(1),w_range(2));  
            end

        end

        // plot divider and fermi level 
        // **scilab use stack, so childern(1) is always the first item

        // plot E_Fermi
        plot(k_range',Ef*ones(tot_k,1),'r:');  
        a=gca();
        a.children(1).children.thickness=1 

        // plot k-divider
        if length(k_div) > 2 then
            ban_max=max(Ek);
            ban_min=min(Ek);
            ban_width=max(Ek)-min(Ek)
            for m=2:length(k_div)-1
                plot(k_range(k_div(m))*ones(10,1),..
                linspace(ban_min-0.05*ban_width,..
                ban_max+0.05*ban_width,10)','k:');
                a.children(1).children.thickness=1
            end
        end

        a=gca();
        // set plot range
        if length(E_bound)==2    
            a.data_bounds=[k_range(1), E_bound(1)...
            ;k_range($), E_bound(2)];
        end
        // set label
        title(desc+' / PBAND-'+string(n),'fontsize',4); 
        ylabel('Energy (eV)',"fontsize", font_size);
        a=gca();
        a.font_size=font_size
        a.tight_limits='on'
        a.thickness=linewidth
        // tick label and number must be tuned simtaneously
        b=a.x_ticks
        b(2)=k_range(k_div);
        b(3)=k_label;
        a.x_ticks=b
    end
endfunction
