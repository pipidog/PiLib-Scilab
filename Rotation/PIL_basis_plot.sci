// **** PURPOSE **** 
// Tplots the rotated coordinate which helps you imagine the 3D story.
// **** VARIABLES ****
// [basis_set]: 3x3, real
// <= classical basis vectors in row.
// [view_point]: string, 'passive' or 'active'
// <= define passive or active rotation
// **** VERSION ****
// 2/23/2014 FIRST BUILT 
// **** COMMENT ****
// Here view_point is meaningless. It's just to show the title.

function []=PIL_basis_plot(basis_set,view_point)
    [lhs,rhs]=argn();
    if rhs==1 then
        view_point='passive';
    end
    if view_point=='passive' then
        basis_set=basis_set';
    end
    scf();set(gca(),"auto_clear","off");
    param3d(linspace(-1,1,2),linspace(0,0,2),linspace(0,0,2))
    e=gce(); e.foreground=color('red'); e.line_style=3; e.thickness=1.0; e.polyline_style=4;;
    param3d(linspace(0,0,2),linspace(-1,1,2),linspace(0,0,2))
    e=gce(); e.foreground=color('gold'); e.line_style=3; e.thickness=1.0; e.polyline_style=4;
    param3d(linspace(0,0,2),linspace(0,0,2),linspace(-1,1,2))
    e=gce(); e.foreground=color('darkgreen'); e.line_style=3; e.thickness=1.0; e.polyline_style=4;

    param3d(linspace(-basis_set(1,1),basis_set(1,1),2),linspace(-basis_set(2,1),basis_set(2,1),2),linspace(-basis_set(3,1),basis_set(3,1),2));
    e=gce(); e.foreground=color('red'); e.thickness=2.0; e.polyline_style=4;
    param3d(linspace(-basis_set(1,2),basis_set(1,2),2),linspace(-basis_set(2,2),basis_set(2,2),2),linspace(-basis_set(3,2),basis_set(3,2),2))
    e=gce(); e.foreground=color('gold'); e.thickness=2.0; e.polyline_style=4;
    param3d(linspace(-basis_set(1,3),basis_set(1,3),2),linspace(-basis_set(2,3),basis_set(2,3),2),linspace(-basis_set(3,3),basis_set(3,3),2))
    e=gce(); e.foreground=color('darkgreen'); e.thickness=2.0; e.polyline_style=4;
    legend(['X','Y','Z']);
    title(view_point);
endfunction
