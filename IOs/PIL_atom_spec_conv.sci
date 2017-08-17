// **** Purpose ****
// This function is to convert the standard atomic species input to 
// code readable row format
// crystal structure. 
// **** Variables ****
// [atom_spec]: nx1, string
// <= the species of atoms. e.g: ['4*Bi','4*O','8*Sr']
// [atom_spec_row]: nx1, string
// => the specie of atom in row: ['Bi,'Bi','Bi','Bi',...]
// **** Version ****
// 09/24/2016: first built
// **** Comment ****
function atom_spec_row=PIL_atom_spec_conv(atom_spec)
    atom_spec_tmp=[];
    for n=1:length(length(atom_spec))
        tmp=strsplit(atom_spec(n));
        tmp_idx=find(tmp=='*');
        tmp_rep=evstr(strcat(tmp(1:tmp_idx-1)));
        tmp_at=strcat(tmp(tmp_idx+1:$));
        atom_spec_tmp=cat(2,atom_spec_tmp,repmat(tmp_at,1,tmp_rep));
    end
    atom_spec_row=atom_spec_tmp;
endfunction

