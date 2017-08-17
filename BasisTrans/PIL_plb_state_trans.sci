// **** Purpose ****
// This function geneates the basis transformation matrix based on the 
// hop.state_info. So PiLab calculation can be handled easier. 
// **** Variables ****
// [state_info]: structure 
// <= hop.state_info variable in PiLab  
// [basis_in]: string, 'c'/'s'/'rc'/'rs'
// <= the input basis
// [basis_out]: string, 'c'/'s'/'rc'/'rs'
// <= the input basis
// [U]: n x n, complex
// => the basis transformation matrix
// **** Version ****
// 11/04/2014
// **** Comment ****
// 1. This function helps you to convert the states in PiLab to other basis.
// However, to make the output consistent, one should make sure you 
// have converted them back after finishing the operations. 
// 2. To use it: A=U_in_out*A*U_in_out'; A=U_in_out*A;
// 3. Remember, to use this, you must make sure you states are complete
// which means you cannot select states in your PiLab_hop input file

function U=PIL_plb_state_trans(state_info,basis_in,basis_out)
    tot_orb=max(state_info(:,3));
    U=[];
    for n=1:tot_orb
        Orb='none';
        L=((length(find(state_info(:,3)==n)))/2-1)/2
        select L
        case 0
            Orb='s'   
        case 1
            Orb='p'    
        case 2
            Orb='d'
        case 3
            Orb='f'
        else
            disp('Error: PIL_plb_state_trans, '...
            +'spin is not included!');
            abort
        end
        [M_out,U_in_out]=PIL_basis_trans(eye(2*(2*L+1),2*(2*L+1))...
        ,Orb,basis_in,basis_out);
        U=PIL_dirsum(U,U_in_out);
    end
endfunction
