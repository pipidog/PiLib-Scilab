// **** Purpose ****
// gnerate a quantum mechanical angular momentun matrix Jx, Jy, Jz
// **** Variables ****
// [J]: 1x1, integer or half-integer
// <= the total angular momentum
// [out_format]: 1x1, 'i' or 'd', default='i'
// <= arrange matrix in increasing or descending order
// [Jx],[Jy],[Jz]: 2J+1 x 2J+1, complex
// => the Jx, Jy, Jz matrix 
// **** Version ****
// 05/01/2014
// **** Comment ****

function [Jx,Jy,Jz]=PIL_J_mat(J,out_format)
    [lhs,rhs]=argn();
    if rhs==1
        out_format='i';
    end
    
    if out_format~='i' & out_format~='d'
        disp('Error! only ''i'' (increasing Mz) or ''d'' (decreasing Mz) are acceptable!');
        pause;
    end  
    
    Jx=zeros(2*J+1,2*J+1);
    Jy=zeros(2*J+1,2*J+1);
    Jz=zeros(2*J+1,2*J+1);
    mz=linspace(-J,J,2*J+1);
    
    for n=1:2*J+1        
        if n~=1 
          Jx(n-1,n)=(1/2)*sqrt((J-mz(n)+1)*(J+mz(n)));
          Jy(n-1,n)=(%i/2)*sqrt((J-mz(n)+1)*(J+mz(n)));
        end
        Jz(n,n)=mz(n);
    end
    Jx=Jx+Jx';
    Jy=Jy+Jy';
     
     // reorder the matrix from +J ~ -J
     
    Jx_tmp=zeros(2*J+1,2*J+1);
    Jy_tmp=zeros(2*J+1,2*J+1);
    Jz_tmp=zeros(2*J+1,2*J+1);

    for n=1:2*J+1
        for m=n:2*J+1
            Jx_tmp(m,n)=Jx(2*J+2-m,2*J+2-n);
            Jy_tmp(m,n)=Jy(2*J+2-m,2*J+2-n);
            Jz_tmp(m,n)=Jz(2*J+2-m,2*J+2-n);
        end
    end
    
    Jx=Jx_tmp+Jx_tmp';
    Jy=Jy_tmp+Jy_tmp';
    Jz=Jz_tmp;
    
    if out_format=='i'
        Jx=PIL_mat_rev(Jx);
        Jy=PIL_mat_rev(Jy);
        Jz=PIL_mat_rev(Jz);
    end         
endfunction
