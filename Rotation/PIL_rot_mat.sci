// **** Purpose ****
// generates the rotation matrix of a classical system
// **** Variables ****
// [Euler_angles]: 1x3, real
// <= the Euler angles
// [view]: 1x1, char: 'passive' or 'active', default='passive'
// <= the view point of your rotation, 
// [conv_seq]: 1x3, char, any 'x','y','z' combintion or 'axis', default='zyz'
// <= define the convention of Euler angle rotation
// [D_CM]: 3x3, real
// => the rotation matrix 
// **** Version ****
// 05/01/2014
// **** Comment ****
// By default, we use 'zyz' convention. Other conventions like 'axes conventions' 
// can be achieved by input 'zyx'. (see 'axes conventions' in wiki.)  
// the input 'view' can be 'active' or 'passive'. In active view, we use
// U(alpha)*U(beta)*U(gamma). In passive view, we use U(gamma)*U(beta)*U(alpha)

// The default values is 'zyz' and 'passive'. In actve, you rotate with a fixed
// basis and in passive view, your rotation axies depends on each roration steps.

// The conv_seq can also input 'axis' to use axis-angle representation. You define
// an axis and its rotation angle. In this case, input the angles in this way:
// [angle_theta,angle_phi,angle_omega], where z=cos(thata), x=sin(thata)*sin(phi), y=sin(theta)*cos(phi)
// and omega is the rotation angle about this unit vector. active=right hand rule, passive=left hand rule

function [D_CM]=PIL_rot_mat(Euler_angles,view,conv_seq)
    [lhs,rhs]=argn();
    select rhs
    case 3
        conv_seq='zyz';
        view='passive';
    case 4
        conv_seq='zyz';
    end

    D_CM=eye(3,3);

    if conv_seq=='axis' then
        rot_axis=[sin(angle_alpha)*cos(angle_beta),sin(angle_alpha)*sin(angle_beta),cos(angle_alpha)];
        cross_mat=[0 -rot_axis(3) rot_axis(2); rot_axis(3) 0 -rot_axis(1); -rot_axis(2) rot_axis(1) 0];
        select view
        case 'active'
            D_CM=eye(3,3)+cross_mat*sin(angle_gamma)+cross_mat^2*(1-cos(angle_gamma));
        case 'passive'
            D_CM=eye(3,3)+cross_mat*sin(-angle_gamma)+cross_mat^2*(1-cos(-angle_gamma));
        end
    else 
        conv_read=msscanf(conv_seq,'%c %c %c');
        for n=1:3 
            select view
            case 'active'
                select conv_read(n)
                case 'x'
                    U=[1 0 0;
                    0 cos(Euler_angles(n)) -sin(Euler_angles(n));
                    0 sin(Euler_angles(n)) cos(Euler_angles(n))];
                case 'y'
                    U=[cos(Euler_angles(n)) 0 sin(Euler_angles(n));
                    0 1 0;
                    -sin(Euler_angles(n)) 0 cos(Euler_angles(n))];
                case 'z'
                    U=[cos(Euler_angles(n)) -sin(Euler_angles(n)) 0;
                    sin(Euler_angles(n)) cos(Euler_angles(n)) 0;
                    0 0 1];
                else
                    disp('Error: PIL_rot_cm');
                    disp('Reason: conv_seq can only have x,y,z!');
                    pause;
                end
                D_CM=D_CM*U; //// alpha, beta, gamma
            case 'passive'
                select conv_read(n)
                case 'x'
                    U=[1 0 0;
                    0 cos(Euler_angles(n)) sin(Euler_angles(n));
                    0 -sin(Euler_angles(n)) cos(Euler_angles(n))];
                case 'y'
                    U=[cos(Euler_angles(n)) 0 -sin(Euler_angles(n));
                    0 1 0;
                    sin(Euler_angles(n)) 0 cos(Euler_angles(n))];
                case 'z'
                    U=[cos(Euler_angles(n)) sin(Euler_angles(n)) 0;
                    -sin(Euler_angles(n)) cos(Euler_angles(n)) 0;
                    0 0 1];
                else
                    disp('Error: PIL_rot_cm');
                    disp('Reason: conv_seq can only have x,y,z!');
                    pause;
                end
                D_CM=U*D_CM; // gamma, beta, alpha
            else
                disp('Error: PIL_rot_cm');
                disp('Reason: view must be active or passive!');
                pause;
            end
        end
    end
endfunction
