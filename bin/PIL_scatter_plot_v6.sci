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
// 08/18/2016 Use "scatter" to adapte to Scilab 6
// **** Comment ****

function PIL_scatter_plot_v6(x,y,c,s,cmap,cmin,cmax)
    [lhs,rhs]=argn();
    select rhs
    case 3
        s=3
        cmap=jetcolormap(64);
        cmin=[];
        cmax=[];
    case 4
        cmap=jetcolormap(64);
        cmin=[];
        cmax=[];
    case 5
        cmin=[];
        cmax=[];
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

    //generate color code
    tot_color=length(cmap(:,1));
    c_div=linspace(cmin,cmax,tot_color+1);
    c_code=zeros(length(x),1);
    for n=1:tot_color
        ind=find(c >= c_div(n) & c < c_div(n+1));
        if ind~=[] then
            c_code(ind,:)=n;
        end
    end
    xset("colormap",cmap);    
    scatter(x,y,s,c_code,"fill");

    a=gca();
    a.tight_limits='on'
    a.font_size=4
    a.thickness=3

    colorbar(cmin,cmax);
    a=gcf();
    a.children(1).font_size=4
    a.children(1).thickness=1
endfunction
//xdel(winsid());
//n=500;
//x=linspace(-%pi,%pi,n)';
//y=sin(x);
//c=linspace(-1,+1,n)';
//s=3;
//cmap=jetcolormap(128);
//PIL_scatter_plot(x,y,c,s,cmap)
