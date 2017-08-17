// **** Purpose ****
// This code performs dampped Fourier transform from t <--> omega. 
// **** Variables ****
// [A]: real, Nx2
// <= input data as function of time. [time, value]. Time must be in 
//    unit of hbar/eV (0.658fs)  
// [E_max]: real, 1x1
// <= max value of FT frequency, in unit of eV
// [dE]: real, 1x1
// <= dE of FT frequency, in unit of eV
// [func]: string, 'sin', 'cos', exp'
// <= def of FT, sin(wt), cos(wt), or exp(iwt)
// [FT_A]: real, Nx2
// => data in frequency domanin. [frequency, value]. omega in unit of eV
// **** Version ****
// 08/10/2015
// **** Comment ****
// This code will output signed FT coefficients. Usually in singal process
// we will only consider its spectrum density, i.e. (abs(FT_A)).^2.

function FT_A=PIL_FT_damp(A,E_max,dE,func)
    [lhs,rhs]=argn();
    if rhs==3 then
        func='sin'
    end
    A=A-repmat(A(1,:),length(A(:,1)),1);
    t=A(:,1);
    dt=t(2)-t(1)

    tot_w=fix(E_max/dE);
    FT_A=zeros(tot_w,2);
    FT_A(:,1)=linspace(0,E_max,tot_w)';
    F_damp=(1-3*(t/t($)).^2+2*(t/t($)).^3)

    for n=1:tot_w    
        select func
        case 'sin'
            FT_A(n,2)=sqrt(1/2*%pi)*sum((sin(FT_A(n,1)*t)).*F_damp.*A(:,2))*dt;
        case 'cos'
            FT_A(n,2)=sqrt(1/2*%pi)*sum((cos(FT_A(n,1)*t)).*F_damp.*A(:,2))*dt;
        case 'exp'
            FT_A(n,2)=sqrt(1/2*%pi)*sum((exp(%i*FT_A(n,1)*t)).*F_damp.*A(:,2))*dt;
        else
            disp('Error: PIL_FT_damp, wrong FT function!')
            abort;
        end
    end
endfunction
