// **** Purpose ****
// transform quantum mechanical basis among atomic orbitals
// **** Variables ****
// [M_in]: n x n, real or complex
// <= input matrix
// [orbital]: 1x1, char: 's', 'p', 'd', 'f'
// <= orbital
// [basis_in]: 1x1, char
// 's'(spherical), 'c'(cubic), 'rc'(LS-sperical), 'rc'(LS-cubic)
// <= input basis
// [basis_out]: 1x1, char
// 's'(spherical), 'c'(cubic), 'rc'(LS-sperical), 'rc'(LS-cubic)
// <= output basis
// [M_out]: nxn, real or complex
// => output matrix after transformation
// [U_in_out]: nxn, real or complex
// => unitary transformatio
// **** Version ****
// 05/01/2014
// **** Comment ****
// 1. Definition of each basis
//  cubic --> J. Phys. C. 13 583
// 2. Notation:
//  U_c_s --> trans a vector form 'c' to 's' representation 
//  ex:U_c_s|c>=|s>
// 3. U_in_out , 
//  M_out=U_in_out*M_in*U_in_out'
//  U_in_out|in> = |out>
// 4. Trick:
//  if you just want the transformation matrix, let M_in=identity matrix
// 5. spinful or spinless
//  this function will automatically detect whether your M_in is spinful or spinless
//  when using 's' or 'c' basis. However, if your M_in is spinless but uses 'rs' or
//  'rc', the code will output error message due to dimensional incosistency.  
// 6. We use PiLib Sparse format to generate the Unitary matrix

function [M_out,U_in_out]=PIL_basis_trans(M_in,orbital,basis_in,basis_out)
    // check input
    if orbital ~= 's' & orbital ~= 'p' & orbital ~= 'd' & orbital ~= 'f' then
        disp('Error: PIL_basis_trans, orbital can be s, p, d, or f !');
        abort
    end
    if basis_in ~= 's' & basis_in ~= 'c' & basis_in ~= 'rs' & basis_in ~= 'rc' then
        disp('Error: PIL_basis_trans, basis_in can be s, c, rs, or rc !');
        abort
    end
    if basis_out ~= 's' & basis_out ~= 'c' & basis_out ~= 'rs' & basis_out ~= 'rc' then
        disp('Error: PIL_basis_trans, basis_out can be s, c, rs, or rc !');
        abort
    end

    // unitary transformation for each oribtal 
    select orbital
    case 's'
        // 's' to 's' (spinless)
        U_s_s=1;

        // 'c' to 's'
        U_c_s=1; 

        // 'rs' to 's'
        U_rs_s=[2,2,0;1,1,1;2,2,1];
        U_rs_s=PIL_sparse(U_rs_s,'sparse','full'); 

        // 'rc' to 'rs'
        U_rc_rs=[2,2,0;1,1,1;2,2,1];
        U_rc_rs=PIL_sparse(U_rc_rs,'sparse','full'); 
        U_rc_s=U_rs_s*U_rc_rs;

    case 'p'
        // 's' to 's' (spinless)
        U_s_s=eye(3,3);

        // 'c' to 's' (spinless) (cheched again on 3/31/2015)
        U_c_s=[3,3,0;1,1,1/sqrt(2);1,2,%i/sqrt(2);3,1,-1/sqrt(2);3,2,%i/sqrt(2); 2 3 1];
        U_c_s=PIL_sparse(U_c_s,'sparse','full'); 

        // 'rs' to 's' (checked again on 3/31/2015)
        U_rs_s=[6,6,0;1,3,1;2,1,1/sqrt(3);2,4,sqrt(2/3);3,2,sqrt(2/3);3,5,sqrt(1/3);4,1,-sqrt(2/3);...
        4,4,sqrt(1/3);5,2,-sqrt(1/3);5,5,sqrt(2/3); 6,6,1];
        U_rs_s=PIL_sparse(U_rs_s,'sparse','full'); 

        // 'rc' to 's'
        U_rc_rs=[6,6,0;1,1,1;2,2,1;3,3,1;4,4,1;5,5,1;6,6,1];
        U_rc_rs=PIL_sparse(U_rc_rs,'sparse','full'); 
        U_rc_s=U_rs_s*U_rc_rs;

    case 'd'
        // 's' to 's' (spinless)
        U_s_s=eye(5,5);

        // 'c' to 's' (spinless)
        U_c_s=[5,5,0;1,1,%i/sqrt(2);1,4,1/sqrt(2);2,2,%i/sqrt(2);2,3,1/sqrt(2);...
        3,5,1;4,2,%i/sqrt(2);4,3,-1/sqrt(2);5,1,-%i/sqrt(2);5,4,1/sqrt(2)];
        U_c_s=PIL_sparse(U_c_s,'sparse','full'); 

        // 'rs' to 's'
        U_rs_s=[10,10,0; 1,5,1; 2,1,sqrt(1/5); 2,6,sqrt(4/5); 3,2,sqrt(2/5); 3,7,sqrt(3/5); 4,3,sqrt(3/5); 4,8,sqrt(2/5);...
        5,4,sqrt(4/5); 5,9,sqrt(1/5); 6,1,-sqrt(4/5); 6,6,sqrt(1/5); 7,2,-sqrt(3/5); 7,7,sqrt(2/5); ...
        8,3,-sqrt(2/5); 8,8,sqrt(3/5); 9,4,-sqrt(1/5); 9,9,sqrt(4/5); 10,10,1];
        U_rs_s=PIL_sparse(U_rs_s,'sparse','full'); 

        // 'rc' to 's'
        U_rc_rs=[10,10,0;1,1,1;2,2,1;3,3,1;4,4,1;5,6,-sqrt(1/6);5,7,-sqrt(5/6);6,5,sqrt(5/6);6,10,sqrt(1/6);7,8,-1;8,9,1;...
        9,6,sqrt(5/6);9,7,-sqrt(1/6);10,5,-sqrt(1/6);10,10,sqrt(5/6)];
        U_rc_rs=PIL_sparse(U_rc_rs,'sparse','full'); 
        U_rc_s=U_rs_s*U_rc_rs;
    case 'f'
        // 's' to 's' (spinless)
        U_s_s=eye(7,7);

        // 'c' to 's'
        U_c_s=[7,7,0;1,2,sqrt(5)/4;1,3,-%i*sqrt(5)/4; 1,5,-sqrt(3)/4; 1,6,-%i*sqrt(3)/4;...
        2,1,%i/sqrt(2); 2,7,1/sqrt(2); 3,2,-sqrt(3)/4;3,3,-%i*sqrt(3)/4;3,5,-sqrt(5)/4;...
        3,6,%i*sqrt(5)/4;4,4,1;5,2,sqrt(3)/4;5,3,-%i*sqrt(3)/4; 5,5,sqrt(5)/4;...
        5,6,%i*sqrt(5)/4; 6,1,-%i/sqrt(2); 6,7,1/sqrt(2); 7,2,-sqrt(5)/4;7,3,-%i*sqrt(5)/4;...
        7,5,sqrt(3)/4;7,6,-%i*sqrt(3)/4 ];
        U_c_s=PIL_sparse(U_c_s,'sparse','full'); 

        //'rs' to 's'
        U_rs_s=[14,14,0;2,1,sqrt(1/7);8,1,-sqrt(6/7);3,2,sqrt(2/7);9,2,-sqrt(5/7);10,3,-sqrt(4/7);4,3,sqrt(3/7);...
        5,4,sqrt(4/7);11,4,-sqrt(3/7);6,5,sqrt(5/7);12,5,-sqrt(2/7);7,6,sqrt(6/7);13,6,-sqrt(1/7);...
        1,7,1;2,8,sqrt(6/7);8,8,sqrt(1/7);3,9,sqrt(5/7);9,9,sqrt(2/7);4,10,sqrt(4/7);10,10,sqrt(3/7);...
        5,11,sqrt(3/7);11,11,sqrt(4/7);6,12,sqrt(2/7);12,12,sqrt(5/7);7,13,sqrt(1/7);13,13,sqrt(6/7);14,14,1];
        U_rs_s=PIL_sparse(U_rs_s,'sparse','full'); 

        // 'rc' to 's'
        U_rc_rs=[14,14,0;2,1,sqrt(5/6);6,1,-sqrt(1/6);1,2,-sqrt(1/6);5,2,sqrt(5/6);1,3,-sqrt(5/6);5,3,-sqrt(1/6);3,4,-1;...
        4,5,1;2,6,sqrt(1/6);6,6,sqrt(5/6);7,7,sqrt(5/12);11,7,sqrt(7/12);10,8,-sqrt(7/12);14,8,-sqrt(5/12);8,9,sqrt(3/4);...
        12,9,-1/2;9,10,1/2;13,10,-sqrt(3/4);7,11,sqrt(7/12);11,11,-sqrt(5/12);9,12,sqrt(3/4);13,12,1/2;...
        8,13,1/2;12,13,sqrt(3/4);10,14,-sqrt(5/12);14,14,sqrt(7/12)];
        U_rc_rs=PIL_sparse(U_rc_rs,'sparse','full'); 
        U_rc_s=U_rs_s*U_rc_rs;
    end

    // transform M_in to 's' basis
    select basis_in
    case 's'
        if prod(size(M_in)==size(U_s_s))==1 then
            U_in_out=U_s_s';
        else
            U_in_out=(PIL_dirsum(U_s_s,U_s_s))';
        end
    case 'c'
        if prod(size(M_in)==size(U_c_s))==1 then
            U_in_out=U_c_s';
        else
            U_in_out=(PIL_dirsum(U_c_s,U_c_s))';
        end
    case 'rs'
        U_in_out=U_rs_s';
    case 'rc'
        U_in_out=U_rc_s';
    end
    
    select basis_out
    case 's'
        if prod(size(M_in)==size(U_s_s))==1 then
            U_in_out=U_in_out*U_s_s;
        else
            U_in_out=U_in_out*(PIL_dirsum(U_s_s,U_s_s));
        end
    case 'c'
        if prod(size(M_in)==size(U_c_s))==1 then
            U_in_out=U_in_out*U_c_s;
        else
            U_in_out=U_in_out*(PIL_dirsum(U_c_s,U_c_s)); 
        end
    case 'rs'
        U_in_out=U_in_out*U_rs_s;
    case 'rc'
        U_in_out=U_in_out*U_rc_s;
    end  
    U_in_out=clean(U_in_out');
    M_out=clean(U_in_out*M_in*U_in_out');
endfunction
