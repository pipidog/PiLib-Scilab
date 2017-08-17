// **** Purpose ****
// gnerate the dy rotation matrix in details
// **** Variables ****
// [J]: 1x1, integer or half-integer
// <= the total angular momentum
// [d_function]: 8 x n, complex
// => the information of each term shown in Saurai p.223
//    [mz1,mz2,r,numerator,denominator,(-1)^(k-m+m')* d_func(4)/d_func(5), 
//     order of cos(B/2), rder of sin(B/2)]
// **** Version ****
// 05/01/2014
// **** Comment ****
// sometimes you may need the dy function of every term in detail. This function
// does it. d-function generator (see Sakurai P.223) 

function [d_function]=PIL_dy_mat(J)
m_z=linspace(J,-J,2*J+1);
count=0;
for m1=1:2*J+1
    for m2=1:2*J+1
        for r=0:2*J+1 // search all allowed summation index, follow Sakurai p.223
            if (J-m_z(m2)-r)>= 0 & (J+m_z(m1)-r) >=0 & (r+m_z(m2)-m_z(m1)) >=0
             count=count+1;
             d_function(count,1)=m_z(m1);
             d_function(count,2)=m_z(m2);
             d_function(count,3)=r; // allowed k
             d_function(count,4)=factorial(J+m_z(m1))*factorial(J-m_z(m1))*factorial(J+m_z(m2))*factorial(J-m_z(m2)); // numerator
             d_function(count,5)=factorial(J-m_z(m2)-r)*factorial(J+m_z(m1)-r)*factorial(r+m_z(m2)-m_z(m1))*factorial(r); // denominator
             d_function(count,6)=(-1)^(r+m_z(m2)-m_z(m1))*(d_function(count,4))^(1/2) / d_function(count,5); // (-1)^(k-m+m')* d_func(4)/d_func(5)
             d_function(count,7)=2*J+m_z(m1)-m_z(m2)-2*r; // order of cos(B/2)
             d_function(count,8)=m_z(m2)-m_z(m1)+2*r; // order of sin(B/2)
            end        
        end
    end
  end
endfunction
