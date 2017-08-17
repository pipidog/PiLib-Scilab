// **** Purpose ****
// calculate the CG coefficients
// **** Variables ****
// [j1_inp],[j2_inp],[j_inp] : 1x1, integer or half-integer
// <= quantum number j1, j2
// [m1_inp],[m2_inp],[m_inp] : 1x1, integer or half-integer
// <= quantum number m1, m2, m
// [cgcoff]: 1x1, real or complex
// => cg coefficient
// **** Version ****
// 05/01/2014
// **** Comment ****
// This program generates the C-G cofficients. 
// Formulas see: "Intro to QM, Tzeng, p238" 

function cgcoff=PIL_cg_coff(j1_inp,j2_inp,m1_inp,m2_inp,j_inp,m_inp)
    j1=j1_inp;  j2=j2_inp;  m1=m1_inp;  m2=m2_inp;  j=j_inp;  m=m_inp
//    if m_inp<0
//        m=-m_inp; m1=-m1_inp; m2=-m2_inp
//    end
//    if j1_inp<j2_inp
//        j1=j2_inp; j2=j1_inp
//    end

    cgcoff=0

    if j_inp<=j1_inp+j2_inp & j_inp>=abs(j1_inp-j2_inp)
        if abs(m1_inp)<=j1_inp & abs(m2_inp)<=j2_inp & abs(m_inp)<=j_inp & (m1_inp+m2_inp)==m_inp
            overall_const=sqrt( (2*j+1)*factorial(j+j1-j2)*factorial(j-j1+j2)*factorial(j1+j2-j)/factorial(j1+j2+j+1))...
            *sqrt(factorial(j+m)*factorial(j-m)*factorial(j1-m1)*factorial(j1+m1)*factorial(j2-m2)*factorial(j2+m2));    
            for k=0:(j1+j2)
                if (j1+j2-j-k)>=0 & (j1-m1-k)>=0 & (j2+m2-k)>=0 & (j-j2+m1+k)>=0 & (j-j1-m2+k)>=0
                    if m==m1+m2
                        cgcoff=cgcoff+overall_const*(-1)^(k)*(factorial(k)*factorial(j1+j2-j-k)*factorial(j1-m1-k)...
                        *factorial(j2+m2-k)*factorial(j-j2+m1+k)*factorial(j-j1-m2+k))^(-1)
                    end
                end
            end
        end
    end



//    if m_inp<0
//        cgcoff=(-1)^(-(j-j1-j2))*cgcoff;
//    end
//    if j1_inp<j2_inp
//        cgcoff=(-1)^(-(j1+j2-j))*cgcoff;
//    end
endfunction
