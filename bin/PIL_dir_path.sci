// **** Purpose ****
//  This function add a '/' to the directary path 
// **** Variables ****
// [dir_path]: 1x1, str
// <= original path 
// => path with '/' ending 
// **** Version ****
// Apr 22, 2016: 1st version
// **** Comment ****

function dir_path=PIL_dir_path(dir_path)
    if dir_path==[] then
        dir_path=pwd();
    end
    tmp=strsplit(dir_path)
    if tmp($)~='/' & tmp($)~='\' then
        dir_path=dir_path+'/';
    end
endfunction
