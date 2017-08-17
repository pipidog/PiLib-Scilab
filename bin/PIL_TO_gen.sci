// **** Purpose ****
// generate a tensor operator basis
// **** Variables ****
// J: 1x1, integer or half-integer
// <= quantu number of total moment
// basis: 1x1, char, 's' or 'c' or 't'
// <= 's'(shperical), 'c'(cubic), 't'(transition)
// T: n x n x n^2, real or complex
// => tensor operator basis
// **** Version ****
// 05/01/2014
// **** Comment ****

function T=PIL_TO_gen(J,basis)
    // check input parameters    
    if J<=0 
        disp('Error! J must be positive!');
    end
    [lhs,rhs]=argn();
    if rhs==1
        basis='s'
    end
    if basis~='s' & basis~='c' & basis~='t'
        disp('Error! only ''s''(shperical), ''c''(cubic), ''t''(transition) are acceptable!');
    end  
    // construct basis
    T=zeros(2*J+1,2*J+1,(2*J+1)^2)

    select basis
    case 't'
        for n=1:(2*J+1)^2
            tmp=zeros((2*J+1)^2,1);
            tmp(n)=1;
            T(:,:,n)=matrix(tmp,2*J+1,2*J+1);
        end
    case 's'
        count=0;
        m_label=linspace(-J,J,2*J+1)
        for K=0:2*J
            for Q=-K:+K
                count=count+1;
                for m1=1:2*J+1
                    for m2=1:2*J+1
                        T(m1,m2,count)=(-1)^(J-m_label(m1))*(2*K+1)^(1/2)*PIL_tri_j_sym(J,J,K,m_label(m2),-m_label(m1),Q);
                    end
                end
            end
        end
    case 'c'
        // spherical
        count=0;
        m_label=linspace(-J,J,2*J+1)
        for K=0:2*J
            for Q=-K:+K
                count=count+1;
                for m1=1:2*J+1
                    for m2=1:2*J+1
                        T(m1,m2,count)=(-1)^(J-m_label(m1))*(2*K+1)^(1/2)*PIL_tri_j_sym(J,J,K,m_label(m2),-m_label(m1),Q);
                    end
                end
            end
        end
        // making KQ label
        O=zeros(T);
        T_label=zeros((2*J+1)^2,2);
        count1=0;
        count2=0;
        for K=0:2*J
            count1=count2+1;
            count2=count2+2*K+1;
            T_label(count1:count2,1)=K;
            T_label(count1:count2,2)=linspace(-K,K,2*K+1)';
        end
        // construct new tensor operators
        count=0;
        for K=0:2*J
            for Q=-K:+K
                count=count+1;
                if Q==0
                    O(:,:,count)=T(:,:,count);
                elseif Q>0
                    pQ_index=find(T_label(:,1)==K & T_label(:,2)==abs(Q));
                    mQ_index=find(T_label(:,1)==K & T_label(:,2)==-abs(Q));
                    O(:,:,count)=(1/sqrt(2))*((-1)^(Q)*T(:,:,pQ_index)+T(:,:,mQ_index));
                elseif Q<0
                    pQ_index=find(T_label(:,1)==K & T_label(:,2)==abs(Q));
                    mQ_index=find(T_label(:,1)==K & T_label(:,2)==-abs(Q));
                    O(:,:,count)=(%i/sqrt(2))*(T(:,:,pQ_index)-(-1)^(Q)*T(:,:,mQ_index));
                end
            end
        end
        T=O;
    end
endfunction
