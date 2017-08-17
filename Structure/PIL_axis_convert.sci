// **** Purpose ****
// convert the crystal axes given in terms of relative angles to 
// Cartesian coordinate representation
// **** Variables ****
// [axis_length]: 1x3, real
// <= the length of a,b,c axis
// [axis_angle]: 1x3, real
// <= the angle between a-b, a-c, b-c
// [output]: 3x3, real
// => coordiante vactros in row
// **** Version ****
// 05/01/2014
// **** Comment ****
// Due to the extra degree of freedom, we will always enforce "c" along [0 0 1] 
// and "a" along [xa,0,za]. Therefore, there are only four variables to solve:xa,
// za, xb,yb; where zb=sqrt(1-xb^2-yb^2);
// a, b, c the length of each axis. 
// Caution: Use PIL_basis_plot(output) to visualize your results as a check ! 
// Note: the parameters axis_length and axis_angle are both 1x3 vectors

function output=PIL_axis_convert(axis_length,axis_angle)
    a=axis_length(:,1); b=axis_length(:,2); c=axis_length(:,3);
    angle_ab=axis_angle(:,1); angle_ac=axis_angle(:,2); angle_bc=axis_angle(:,3);  
    Zc=c;
    Za=a*c*cos(angle_ac)/Zc;
    Xa=sqrt(a^2-Za^2);
    Zb=b*c*cos(angle_bc)/Zc;
    Xb=(a*b*cos(angle_ab)-Za*Zb)/Xa;
    Yb=sqrt(b^2-Xb^2-Zb^2);
    output=[Xa 0 Za;Xb Yb Zb; 0 0 Zc];
endfunction
