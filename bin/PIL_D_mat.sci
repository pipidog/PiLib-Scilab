// **** Purpose ****
// gnerate a quantum mechanical rotation matrix
// **** Variables ****
// [J]: 1x1, integer or half-integer
// <= the total angular momentum
// [Euler_angles]: 1x3, real
// <= Euler angles
// [basis]: 1x1, char, 's'(shperical)or 'c'(cubic), default='s'
// <= the output basis
// [out_format]: 1x1, 'i' or 'd', default='i'
// <= arrange matrix in increasing or descending order
// [D_mat]: 2J+1 x 2J+1, complex
// => the rotation matrix 
// **** Version ****
// 05/01/2014
// **** Comment ****
// D-function generator (see Tzeng V-II p353 ). The input argunment 'basis'
// can be 's' (shperical harmonics) or 'c' (cubic harmonics). 

function [D_matrix]=PIL_D_mat(J,Euler_angles,basis,out_format)
    angle_alpha=Euler_angles(1);
    angle_beta=Euler_angles(2);
    angle_gamma=Euler_angles(3);
    [lhs,rhs]=argn();
    select rhs
    case 2
        basis='s';
        out_format='i';
    case 3
        out_format='i';
    end

    if basis~='s' & basis~='c'
        disp('Error! only ''s'' (shperical harmonics) or ''c'' (tesseral harmonics) are acceptable!');
        pause;
    end  

    if out_format~='i' & out_format~='d'
        disp('Error! only ''i'' (increasing Mz) or ''d'' (decreasing Mz) are acceptable!');
        pause;
    end  

    m_label=linspace(-J,J,2*J+1);
    D_matrix=zeros(2*J+1,2*J+1);
    D_y=zeros(2*J+1,2*J+1);
    for m1=1:2*J+1
        for m2=1:2*J+1
            // search all allowed summation index, follow Tzeng V-II Ed.3
            term1=sqrt(factorial(J+m_label(m1))*factorial(J-m_label(m1))*factorial(J+m_label(m2))*factorial(J-m_label(m2)));
            for r=0:2*J+1 
                if (J-m_label(m2)-r)>= 0 & (J+m_label(m1)-r) >=0 & (r+m_label(m2)-m_label(m1)) >=0
                    term2=((-1)^(r)*factorial(J-m_label(m2)-r)*factorial(J+m_label(m1)-r)*factorial(r+m_label(m2)-m_label(m1))*factorial(r))^(-1); 
                    term3=((cos(angle_beta/2))^(2*J+m_label(m1)-m_label(m2)-2*r))*((-sin(angle_beta/2))^(m_label(m2)-m_label(m1)+2*r));
                    D_y(m2,m1)=D_y(m2,m1)+term1*term2*term3; 
                end        
            end
            D_matrix(m2,m1)=exp(-%i*m_label(m2)*angle_alpha)*D_y(m2,m1)*exp(-%i*m_label(m1)*angle_gamma); 
        end
    end

    //change basis
    if basis=='s'
        D_matrix=clean(D_matrix);
    elseif basis=='c'
        U=eye(2*J+1,2*J+1);
        for n=1:2*J+1
            m=m_label(n);
            index_p=find(m_label==abs(m));
            index_m=find(m_label==-abs(m));
            if m>0
                U(index_p,n)=(1/sqrt(2))*(-1)^(abs(m));
                U(index_m,n)=(1/sqrt(2));
            elseif m<0
                U(index_p,n)=(%i/sqrt(2));
                U(index_m,n)=-(%i/sqrt(2))*(-1)^(abs(m));
            end
        end   
        D_matrix=clean(U'*D_matrix*U); 
    end

    //change output format
    if out_format=='i'
        D_matrix=D_matrix;
    elseif out_format=='d'
        D_matrix=PIL_mat_rev(D_matrix);
    end 
endfunction

// confirm D in sh & D in oh get the same result:
//
//Q=PIL_tensor_op(1,'o');
//T=PIL_tensor_op(1,'s');
//DO=PIL_D_mat(2,0,%pi/3,0,'o');
//DT=PIL_D_mat(2,0,%pi/3,0,'s');
//tmp1=DO*[0 0 0 0 1]';
//tmp2=DT*(1/sqrt(2))*((-1)^(2)*[0 0 0 0 1]'+[1 0 0 0 0]');
//T1=zeros(3,3);
//T2=zeros(3,3);
//for n=1:5
//    T1=T1+tmp1(n)*O(:,:,n+4);
//    T2=T2+tmp2(n)*T(:,:,n+4);
//end
//// T1=T2 !
