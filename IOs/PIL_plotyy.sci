// **** Purpose ****
// generate a double y axes plot
// **** Variables ****
// [x]: 1 x n, real
// <= x data
// [y1],[y2]: 1 x n, real
// <= y1,y2 data
// [style1,sytle2]: string
// <= the same style input as the system default
// **** Version ****
// 05/01/2014
// **** Comment ****

function PIL_plotyy(x,y1,y2,style1,style2)
    [lhr,rhs]=argn();
    select rhs
    case 3
        style1='-';
        style2='-';
    case 4
        style2='-';
    end
    clf();
    subplot(211);
    plot(x,y1,style1);  // your first figure
    a1 = gca();
    subplot(212)
    plot(x,y2,style2);  // your second figure
    a2 = gca();
    a2.axes_visible = ["off", "on","on"];
    a2.y_location ="right";
    a1.axes_bounds=[0 0 1 1];  // modify the first figure to occupy the whole area
    a2.axes_bounds=[0 0 1 1]; // modify the second figure to occupy the whole area
    a2.filled = "off";  
endfunction

