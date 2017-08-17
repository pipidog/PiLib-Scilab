// **** Purpose ****
// generate the spin-orbital coupling
// **** Variables ****
// [orbital]: (3x1, real) or (char)
// <= (j,l,s) or ('s'/'p'/'d'/'f')
// [coupling]: 1x1, real 
// <= the coupling constant
// [basis_out]: 1x1, char: 's', 'c', 'rs', 'rc'
// <= if orbital is ('s'/'p'/'d'/'f'), you can specify it output basis
// [LS_Mat]: N x N, real
// => spin-orbital coupling matrix
// **** Version ****
// 05/01/2014
// **** Comment ****
// 1. This function allows you generatess the LS coupling by manually input the
//    (j,l,s) quantum number or ('s'/'p'/'d'/'f') orbital. 
// 2. if you use ('s'/'p'/'d'/'f') orbital, you can specify its output basis

function LS_Mat=PIL_LS_coup(orbital,coupling,basis_out)
    // give a (j,l,s), it generates the unit LS matrix in 'rs' basis
    function LS_Mat_jls=LS_coup_jls(j,l,s)
        if (j>abs(l+s)) | (j<abs(l-s)) then
            disp('Error: PIL_LS_coup, j is not allowed!');
            abort
        end
       LS_Mat_jls=(j*(j+1)-l*(l+1)-s*(s+1))/2*eye(2*j+1,2*j+1);
    endfunction
    
    [lhs,rhs]=argn();
    select length(orbital)
    case 3
        LS_Mat=LS_coup_jls(orbital(1),orbital(2),orbital(3));
    case 1

        if rhs==2 then
            basis_out='rs';
        end

        select orbital
        case 's'
            LS_Mat=LS_coup_jls(1/2,0,1/2);
        case 'p'
            LS_Mat=PIL_dirsum(LS_coup_jls(1/2,1,1/2),LS_coup_jls(3/2,1,1/2));
        case 'd'
            LS_Mat=PIL_dirsum(LS_coup_jls(3/2,2,1/2),LS_coup_jls(5/2,2,1/2));
        case 'f'
            LS_Mat=PIL_dirsum(LS_coup_jls(5/2,3,1/2),LS_coup_jls(7/2,3,1/2));
        end
        [LS_Mat,U_in_out]=PIL_basis_trans(LS_Mat,orbital,'rs',basis_out);
    end
    LS_Mat=coupling*LS_Mat
endfunction
