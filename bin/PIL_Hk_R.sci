// **** Purpose ****
// This function generates H(k) from hr.dat of wannier 90. 
// **** Variables ****
// [k]: 1x3, real
// <= [kx,ky,kz]
// [lat_vec]: 3x3, real
// <= lattice row vectors 
// [uc_index]: n x 3 , int
// <= unit cell index, n1*a1+n2*a2+n3*a3
// [uc_deg]:  tot_uc_index x 1, int / []
// <= degenerancy of each uc_index, if [], then assume 1 for all
// [Hr]: tot_state x tot_state x tot_uc_index, complex
//       tot_state x tot_state*tot_uc_index, complex or complex-sparse
// <= <n(0)|Hr|m(R)> matrix, can be 2 or 3 dimensional form. 
// [Hk]: tot_state x tot_state, complex 
// => H(k)
// [err]: 1x1, int
// => if H(k) doesn't pass hermitian check, err=1, otherwise, err=0. 
// **** Version ****
// 01/16/2016 first built. fully rewrite to use new notations.
// 03/25/2016 Hr can be in 2 dimensional form to speed up the code. 
//            (high dimension matrix is much slower!)
// **** Comment ****


function [Hk,err]=PIL_Hk_R(k,lat_vec,uc_index,uc_deg,Hr)
    if uc_deg==[] then
        uc_deg=ones(uc_index(:,1));
    end    

    // check Hr format
    tot_uc=length(uc_index(:,1));
    if ndims(Hr)==3 then
        tot_state=length(Hr(:,1,1));
        Hr=matrix(Hr,tot_state,tot_state*tot_uc);
    end
    tot_state=length(Hr(:,1));
    
    // construct Hk
    Hk=zeros(tot_state,tot_state);
    for n=1:tot_uc
        if max(abs(Hr(:,(n-1)*tot_state+1:n*tot_state)))>=1e-5 then
            R=uc_index(n,:)*lat_vec;
            Hk=Hk+Hr(:,(n-1)*tot_state+1:n*tot_state)*exp(%i*k*R')/uc_deg(n);
        end
    end

    // check hermitian
    hermitian_err=max(abs(Hk-Hk'));
    if hermitian_err<1e-4 then
        Hk=(Hk+Hk')/2;
        err=0
    else
        disp('Error: PIL_Hk_R, Hk is not Hermitian!')
        disp(cat(2,'hermitian_err=',string(hermitian_err)));
        disp(cat(2,'max abs value of Hk=',string(max(abs(Hk)))));
        err=1
        Hk=(Hk+Hk')/2
    end
endfunction
