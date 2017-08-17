// **** Purpose ****
// This code generates the primitive cell and the sublattices of a super
// cell, i.e, enlarge the unit cell of a given system. 
// **** Variables ****
// [primitive]: 3x3, real
// <= primitive vectors in row vectors
// [sublatt_format]: 1x1, string, 'coordinate' / 'coefficient'
// <= specify the format of the sublattice position, coefficients of
//    primitive vectors or cartesian 
// [supercell]: 1x3, integer
// <= the size of the supercell along each direction
// [sup_primitive]: 3x3, real
// [screen_output]: 1x1, string, 'on' / 'off'
// <= whether to display formatted results in console, so one can use it
//    to ab initio inputs directly.  
// => primitive cell of the supercell
// [sup_sublatt]: Nx3, real
// => sublattice of the supercell
// **** Version ****
// 04/21/2015 1st version
// **** Comment ****
// 1. This function is useful to generate supercell coordinate especially
// for ab initio calculations.
// 2. Note: this fucntion assumes all coordinates are in 3 dimention

function [sup_primitive,sup_sublatt]=PIL_supercell(primitive,sublatt...
    ,supercell,sublatt_format,screen_output)
    [lhs,rhs]=argn();
    select rhs
    case 3
        sublatt_format='coordinate'
        screen_output='off'
    case 4
        screen_output='off'
    end

    // define parameters ================================================
    tot_unit_sublatt=length(sublatt(:,1));
    tot_cell=prod(supercell);
    tot_sublatt=tot_unit_sublatt*tot_cell;
    // generate all sublattice coordinate ===============================
    sup_sublatt=zeros(tot_sublatt,3);
    count=0;
    for n1=1:supercell(1)
        for n2=1:supercell(2)
            for n3=1:supercell(3)
                select sublatt_format
                case 'coordinate'
                    new_cell_pos=(n1-1)*primitive(1,:)..
                    +(n2-1)*primitive(2,:)+(n3-1)*primitive(3,:);
                case 'coefficient'
                    new_cell_pos=[(n1-1),(n2-1),(n3-1)];
                end
                
                count=count+1;
                r1=(count-1)*tot_unit_sublatt+1;
                r2=r1+tot_unit_sublatt-1;
                sup_sublatt(r1:r2,:)=repmat(new_cell_pos...
                ,tot_unit_sublatt,1)+sublatt;
            end
        end
    end
    // output new primitive cell
    sup_primitive=zeros(primitive);
    for n=1:3
        sup_primitive(n,:)=supercell(n)*primitive(n,:);
    end
    
    select screen_output
    case 'on'
        disp('super primitive cell');
        mfprintf(6,'%7.4f %7.4f %7.4f\n',sup_primitive);
        disp('super sublattice (total='+string(tot_sublatt)+')')
        mfprintf(6,'%7.4f %7.4f %7.4f\n',sup_sublatt);
    end
endfunction

////example: graphene supercell      
//primitive=[1/2,sqrt(3)/2,0; -1/2, sqrt(3)/2,0; 0, 0, 10];            
//sublatt=[0,0,0;0,-sqrt(3)/3,0]
//supercell=[7,7,1];
//[sup_primitive,sup_sublatt]=PIL_supercell(primitive,sublatt,supercell)

