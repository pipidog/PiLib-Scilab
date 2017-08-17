// **** Purpose ****
// plot (x,y) with z-value labeled by color
// **** Variables ****
// [x]: nx1, real
// x-value
// [y]: nx1, real
// y-value
// [c]: nx1, real
// weighting value of each (xi,yi)
// [s]: 1x1, int
// size of the scatter
// [cmap]: nx3, real
// colormap, e.g. oceancolormap(64)
// <= description of the variable 
// [cmax]: 1x1, real
// <= if you want the color to use the same min and max to plot, set cmax
//    then color will be divided 
// **** Version ****
// Apr 17, 2016 first built
// 06/03/2016 add cmax parameter
// **** Comment ****
function PIL_scatter_plot(x,y,c,s,cmap,cmin,cmax)
    [lhs,rhs]=argn();
    select rhs
    case 3
        s=5
        cmap=[];
        cmin=[];
        cmax=[];
    case 4
        cmap=[];
        cmin=[];
        cmax=[];
    case 5
        cmin=[];
        cmax=[];
    end
    
    if cmap==[] then
        c_map_low=[0.75 0.75 0.75]
        c_map_high=[0.0 0.0 0.6]
        cmap=flipdim(PIL_k_path([c_map_high;c_map_low],64),1);
    end
    if cmin==[] then
        cmin=min(c);
    end
    if cmax==[] then
        cmax=max(c);
    end
    // check variables
    if cmin > min(c) then
        disp('Error: PIL_scatter_plot, cmin < min(c)');
        abort
    end
    if cmax < min(c) then
        disp('Error: PIL_scatter_plot, cmax < max(c)');
        abort
    end

    xset("colormap",cmap);
    tot_color=length(cmap(:,1));
    c_div=linspace(cmin,cmax,tot_color+1);
    for n=1:tot_color
        ind=find(c >= c_div(n) & c < c_div(n+1));
        if ind~=[] then
            plot(x(ind),y(ind),'.','MarkerSize',s,'color',cmap(n,:));
            a=gce();
            a.children.line_mode='off'
        end
    end
    colorbar(cmin,cmax);
    a=gca();
    a.tight_limits='on'
    a.font_size=3
    a.thickness=3

    a=gcf();
    a.children(1).font_size=3
    a.children(1).thickness=1
    a.background=-2;
endfunction
//
//clear; xdel(winsid());
//n=500;
//x=linspace(-%pi,%pi,n)';
//y=sin(x);
//c=linspace(-1,+1,n)';
//s=3;
//cmap=jetcolormap(128);
//PIL_scatter_plot(x,y,c,s,cmap)
