// **** Purpose ****
// find Euler angles between two coordinates
// **** Variables ****
// [basis_set_1],[basis_set_2]: 3x3, real
// <= the basis vectors in row
// [view]: 1x1, string, 'passive' or 'active'
// <= the view point of your rotation
// [conv_seq]: string, 'x,y,z' combination or 'axis'
// <= the convention of the Euler angles
// [Euler_angles]: 1x3, real
// => the Euler angles
// **** Version ****
// 05/01/2014
// **** Comment ****
// Euler angles finder between [x1,y1,z1] and [x2,y2,z2].
// A reference can be found by seraching keyword:"Rotation formalisms in three dimensions" in wiki.
// However, in their foumula, they still cannot elimite the arbitariness Euler angles. 
// So we calculate all possible choice and make sure they are identical.

// Note: 'basis_set_1' and 'basis_set_2' are both 3x3 matrices. The should construct by [x,y,z] and 
// [x',y',z']. If you just want to the the Eulers angles of a particular rotation matrix, set:
// basis_set_1=eye(3,3) and basis_set_2=D.

//Here the parameters 'view' and and 'conv_seq' follows the same definition as PIL_rot_cm 

function [Euler_angles]=PIL_Euler_finder(basis_set_1,basis_set_2,view,conv_seq)
    [lhs,rhs]=argn()
    if rhs==2 then
        view='passive';
    end
    if view~='passive' & view~='active' then
        disp('Error: PIL_Euler_finder');
        disp('Reason: view can only be passive or active!');
    end
    if sum(size(basis_set_1)==[3,3])~=2 | sum(size(basis_set_2)==[3,3])~=2 then
        disp('Errir! basis vectors must be 3x1 column vectors!');
    end
    
    // generate U matrix ============================
    for n=1:3
        for m=1:3
            U(n,m)=basis_set_1(:,n)'*basis_set_2(:,m);
        end
    end
    // coaxis =======================================
    N=U(3,2)*basis_set_2(:,1)-U(3,1)*basis_set_2(:,2);
    N=N/norm(N);
    // list all possible Euler angles ===============
    angle_alpha=[acos(N'*basis_set_1(:,2)),acos(-N'*basis_set_1(:,2)),-acos(N'*basis_set_1(:,2)),-acos(-N'*basis_set_1(:,2))];
    angle_beta=acos(basis_set_2(:,3)'*basis_set_1(:,3));
    angle_gamma=[acos(N'*basis_set_2(:,2)),acos(-N'*basis_set_2(:,2)),-acos(N'*basis_set_2(:,2)),-acos(-N'*basis_set_2(:,2))];

    count=0;
    for n=1:4
        for m=1:4
            Euler_angles_tmp=[angle_alpha(n),angle_beta,angle_gamma(m)];
            D=PIL_rot_cm(Euler_angles_tmp(1),Euler_angles_tmp(2),Euler_angles_tmp(3),view,conv_seq);
            if sum(abs(clean(D-U)))==0 then
                count=count+1;
                Euler_angles(count,:)=Euler_angles_tmp;
            end
        end
    end
    if sum(abs(clean(basis_set_1-basis_set_2)))~=0 & sum(abs(Euler_angles))==0 then
        disp('Error: PIL_Euler_finder');
        disp('Reason: cannot find Euler angles! Try to change view!');
    end
endfunction
