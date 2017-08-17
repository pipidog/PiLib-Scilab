// **** Purpose ****
// This function can generate an xsf file for visualization of given
// crystal structure. 
// **** Variables ****
// [file_title]: 1x1, string
// <= the file name without extension. extension will alwaus be xsf.
// [lat_vec]: 3x3, real
// <= the lattice row vectors
// [sublat_label]: nx1 or 1xn, integer / char
// <= the atom type of each sublattice. can be atomic number or atom name
//    e.g. 73 or 'Ta'
// [sublat_pos]: nx3, real
// <= the cartisian x,y,z positions of each sublattice. 
// **** Version ****
// 02/01/2016: first built
// **** Comment ****

function PIL_crystal_xsf(file_title,lat_vec,sublat_label,sublat_pos)
    tmp=size(sublat_label);
    if tmp(2)~=1 then
        sublat_label=sublat_label';
    end
    
    fid=mopen(file_title+'.xsf','w')
    mfprintf(fid,'%s\n','CRYSTAL')
    mfprintf(fid,'  %s\n','PRIMVEC')
    mfprintf(fid,'    %9.5f  %9.5f  %9.5f\n',lat_vec)
    mfprintf(fid,'  %s\n','PRIMCOORD')
    mfprintf(fid,'    %d  1\n',length(sublat_pos(:,1)))
    for n=1:length(sublat_pos(:,1))
        select type(sublat_label)
        case 1
            mfprintf(fid,'    %3d  %9.5f  %9.5f  %9.5f\n',..
            sublat_label(n),sublat_pos(n,:));
        case 10
            mfprintf(fid,'    %2s  %9.5f  %9.5f  %9.5f\n',..
            sublat_label(n),sublat_pos(n,:));
        end
    end
    mclose(fid)
endfunction
