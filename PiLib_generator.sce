clc; clear; 
// Main ================================================================
bin_folder=PiLib_path+'bin'
folder_name=['BasisTrans','Hamiltonian','PiLab','IOs','Manybody','Math',..
'Property','Rotation','Structure','TensorOp','Utility'];

tic()
for n=1:length(length(folder_name))
    f_list=dir(PiLib_path+folder_name(n));
    for m=1:length(length(f_list(2)))
        if grep(f_list(2)(m),'.sci')~=[] & grep(f_list(2)(m),'.sci~')==[]
            copyfile(PiLib_path+folder_name(n)+'/'+f_list(2)(m),bin_folder);
        end
    end
end
exec(PiLib);
