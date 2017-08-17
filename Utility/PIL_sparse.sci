// **** PURPOSE **** 
// This function generates the sparse matrix from a full matrix or vice versa 
// **** VARIABLES ****
// [A]: (mx3) / (nxn)
// <= the matrix, either full or sparse
// [inp_type]: string, 'sparse' or 'full'
// <= specify your input
// [val_filter]: 1x1, real, default=10^-6
// <= input a value to fliter small values
// [inp_zone]: string, 'tri' or 'all'
// <= if 'tri', sparse matrix only records the upper triangle
//    if 'all', sparse matrix only records all non-zeros  
//    if hermitian matrix, use 'tri' can further save your storage.
// [B]: (mx3) / (nxn)
// => output sparse or full matrix 
// **** VERSION ****
// 2/23/2014 FIRST BUILT 
// 5/17/2014 fix problem when complex values
// 6/4/2014 fix bug, when indices are not integer type
// 4/26/2014 fix bug when input is not a square matrix
// **** COMMENT ****

function [B]=PIL_sparse(A,inp_type,inp_zone,val_fliter)
    [lhs,rhs]=argn();
    select rhs
    case 2 
        inp_zone='all'
        val_fliter=1e-7;
    case 3
        val_fliter=1e-7;
    end

    select inp_type
    case 'sparse'
        B=zeros(real(A(1,1)),real(A(1,2)));
        for n=2:length(A(:,1))
            B(round(real(A(n,1))),round(real(A(n,2))))=A(n,3);
        end
        select inp_zone
        case 'tri'
            B=B+B'-diag(diag(B))
        case 'all'

        end
    case 'full'
        A_dim=size(A);        
        B=[A_dim(1),A_dim(2),0];
        for n=1:A_dim(1)
            A_nz=find(abs(A(n,:))>=val_fliter);
            for m=1:length(A_nz)
                select inp_zone
                case 'tri'
                    if A_nz(m) >=n then
                        B=cat(1,B,[n,A_nz(m),A(n,A_nz(m))]);
                    end
                case 'all'
                    B=cat(1,B,[n,A_nz(m),A(n,A_nz(m))]);
                else
                    disp('Error: PIL_sparse, inp_zone can be ''tri'' or ''all'' only!');
                    abort  
                end
            end
        end
    else 
        disp('Error: PIL_sparse, inp_type can be ''full'' or ''sparse'' only!')
        abort
    end
endfunction
