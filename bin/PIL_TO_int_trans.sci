// **** Purpose ****
// calculate a multipolar exchange matrix to different super basis
// **** Variables ****
// A_in: n1^2 x n2^2, real or complex
// <= The input exchange coupling matrix
// TO1_in: (n1 x n1 x n1^2, real or comple ) / (char: 's', 'c', 't')
// <= The input super basis of site 1
// TO2_in: (n2 x n2 x n2^2, real or comple ) / (char: 's', 'c', 't')
// <= The input super basis of site 2
// TO1_out: (n1 x n1 x n1^2, real or comple ) / (char: 's', 'c', 't')
// <= The output super basis of site 1 
// TO2_out: (n2 x n2 x n2^2, real or comple ) / (char: 's', 'c', 't')
// <= The output super basis of site 2  
// A_out: n1^2 x n2^2, real or complex
// => The output exchange coupling matrix
// **** Version ****
// 05/01/2014
// **** Comment ****
// This code transforms a given interaction defined in one tensor operator basis
// to another tensor operator basis. Since the angular momentums can be different
// on each s. It is possible that two ss have different operator basis. So,
// in this function, you need to input all the basis of each s. basis1 means
// the operators basis represents your row index in the A and TO1_ for 
// column index. Those basis parematers can be 'c','s' or given Tensor operator 
// basis. A_out=U1*A_in*U2.'

function [A_out,U1,U2]=PIL_TO_int_trans(A_in,TO1_in,TO2_in,TO1_out,TO2_out)
    s1_basis_op_num=length(A_in(:,1));
    s2_basis_op_num=length(A_in(1,:));

    select type(TO1_in)
    case 10
        T1_in=PIL_TO_gen((sqrt(s1_basis_op_num)-1)/2,TO1_in);
    case 17
        T1_in=TO1_in;
    end

    select type(TO2_in)
    case 10
        T2_in=PIL_TO_gen((sqrt(s2_basis_op_num)-1)/2,TO2_in);
    case 17
        T2_in=TO2_in;
    end

    select type(TO1_out)
    case 10
        T1_out=PIL_TO_gen((sqrt(s1_basis_op_num)-1)/2,TO1_out);
    case 17
        T1_out=TO1_out;
    end

    select type(TO2_out)
    case 10
        T2_out=PIL_TO_gen((sqrt(s2_basis_op_num)-1)/2,TO2_out);
    case 17
        T2_out=TO2_out;
    end

    U1=zeros(s1_basis_op_num,s1_basis_op_num);
    U2=zeros(s2_basis_op_num,s2_basis_op_num);
    for n=1:s1_basis_op_num
        for m=1:s1_basis_op_num
            U1(n,m)=sum(diag(T1_out(:,:,n)'*T1_in(:,:,m)));
        end
    end
    for n=1:s2_basis_op_num
        for m=1:s2_basis_op_num
            U2(n,m)=sum(diag(T2_out(:,:,n)'*T2_in(:,:,m)));
        end
    end

    U1=clean(U1);
    U2=clean(U2);
    A_out=clean(U1*A_in*U2.');

endfunction

//example:
//A_in=eye(9,9);
//A_out=mylib_MH_int_trans(A_in,'t','t','h','h');
// or:
//T1=PIL_TO_gen(1,'t');
//T2=PIL_TO_gen(1,'h');
//Vout=mylib_MH_int_trans(A_in,T1,T1,T2,T2);
