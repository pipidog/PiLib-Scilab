// **** Purpose ****
// This function geneates the basis transformation matrix based on the 
// flq.state_info. So PiLab calculattion can be handled easier. 
// **** Variables ****
// [state_info]: structure 
// <= flq.state_info variable in PiLab  
// [basis_in]: string, 'c'/'s'/'rc'/'rs'
// <= the input basis
// [basis_out]: string, 'c'/'s'/'rc'/'rs'
// <= the input basis
// [U]: n x n, complex
// => the basis transformation matrix
// **** Version ****
// 11/04/2014
// **** Comment ****
// 1. This function helps you to convert the floquet states in PiLab to 
// other basis. However, to make the output consistent, one should make 
// sure you have converted them back after finishing the operations. 
// 2. This function convert flq.state_info into hop.state_info format.
// Then call PIL_plb_state_trans and direct sum to correct block.
// 3. To use it: A=U_in_out*A*U_in_out'; A=U_in_out*A;

function U=PIL_plb_flq_state_trans(state_info,basis_in,basis_out)
    flq_order=max(state_info(:,2));
    tot_state=round(length(state_info(:,1))/(2*flq_order+1));
    state_info=state_info(1:tot_state,[1,3:6]);

    U=PIL_plb_state_trans(state_info,basis_in,basis_out)
    U_in_out=[];
    for n=1:2*flq_order+1
        U_in_out=PIL_dirsum(U_in_out,U);
    end
endfunction

//clear; exec(PiLib);
//cd('C:\Users\pipidog\Dropbox\My Project\PiLab projects\Graphene_sp');
//PiLab_loader('Graphene_sp','flq');
//load('Graphene_sp_flq.sod');
//state_info=flq.state_info;
//basis_in='c';
//basis_out='rc';
//
