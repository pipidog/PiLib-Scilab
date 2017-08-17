// **** Purpose ****
// This code calculate the Floquet Bessel integrals for a homogeneous 
// (no direction) perodic electric field. 
// **** Variables ****
// [w]: 1x1, real
// <= frequency
// [p]: 1x1, integer
// <= photon process order
// [A_amp]: 1x1 / 1x2 / 1x3 , real
// <= vector potential amplitude Ax, Ay, Az
// [A_phase]: 1x1 / 1x2 / 1x3, real
// <= Phase term of Ax, Ay, Az
// [r]: 1x1 / 1x2 /1x3, real
// <= position vector to perform integration
// **** Version ****
// 7/9/2014 First Built
// **** Comment ****
// 1. About the theory
// This code requires the vector potential independent of space. For 
// 1D and 2D, it is easy to achieve by considering a electric field 
// along z-direction. So, for 2D or 1D, we can assume z=0. For 3D case,
// it can not be easy achieved (maybe long wave length) 
// For calculation details, see PRL 110 200403
// For the form of E field, see PRL 112 156801
// and also my note for Floquet theory
// 2 about the input
// This code accepts 1,2,3 dimension inputs (A_amplitude, A_phase, r)
// but eventually, it will automatically transform to 3D for calculation.
// 3. About this function
// The form of vector potential
// A=[Ax*sin(w*t+s1),Ay*sin(wt+s2),Az*sin(wt+s3)], r=[rx,ry,rz]
// Floquet integral:
// J(w,p,[Ax,Ay],phi,r)=(1/T)*integral( exp(%i*(A*r-pwt) )dt, t=0~T
// 4. For 2D case (Az=0) with Ax=Ay, then
// (s2-s1)=0 --> linear (right up to left down)
// (s2-s1)=%pi/2 --> circular
// (s2-s1)=%pi --> linear (left up to right down)   

function F_val=PIL_flqint(w,p,A_amp,A_phase,r)
    A_amp=PIL_vec_3d(A_amp);
    A_phase=PIL_vec_3d(A_phase);
    r=PIL_vec_3d(r);
    
    T=2*%pi/w;
    func=string(A_amp(1)*r(1))+'*sin('+string(w)+'*t+'+string(A_phase(1))+')+'...
        +string(A_amp(2)*r(2))+'*sin('+string(w)+'*t+'+string(A_phase(2))+')+'...
        +string(A_amp(3)*r(3))+'*sin('+string(w)+'*t+'+string(A_phase(3))+')';
    func_form='y=(1/'+string(T)+')*exp(%i* ('+func+'-'+string(p*w)+'*t)  )';
    deff('y=f(t)',func_form);
    F_val=clean(intc(0,T,f));
endfunction
