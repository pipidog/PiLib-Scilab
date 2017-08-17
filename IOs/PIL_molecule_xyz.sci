// **** Purpose ****
// This function output the input coordinate to xyz format
// **** Variables ****
// [file_title]: 1x1, string
// <= the file name without extension. extension will always be xyz. 
// [atom_type]: 1xn / nx1, string or number
// <= atom type of each coordinate
// [atom_pos]: nx3, real
// <= coordinate of each atom. In cartisian coordinate 
// **** Version ****
// 04/09/2016 first built
// **** Comment ****
function PIL_molecule_xyz(file_title,atom_type,atom_pos)
    tmp=size(atom_type);
    if tmp(2)~=1 then
        atom_type=atom_type';
    end
    
    fid=mopen(file_title+'.xyz','w');
    mfprintf(fid,'%d\n',length(atom_pos(:,1)));
    mfprintf(fid,'molecule structure\n');
    mfprintf(fid,'%2s   %11.7f   %11.7f   %11.7f\n'..
    ,atom_type,atom_pos)
    mclose(fid);
endfunction
