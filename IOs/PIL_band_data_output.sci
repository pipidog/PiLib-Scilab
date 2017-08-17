// This function is to convert band data to format that is easy for 
// origin plotting. 
// Ek_in must be arranged in tot_band x tot_k x tot_proj format.

function PIL_band_data_output(file_name,Ek_in)
    Ek_size=size(Ek_in);
    if length(Ek_size)==2 then
        Ek_size=[Ek_size,1];
    end
    tot_ban=Ek_size(1); 
    tot_k=Ek_size(2);
    tot_proj=Ek_size(3);
    
    Ek_out=zeros(tot_ban*tot_k,tot_proj);
    for n=1:tot_proj
        Ek_tmp=Ek_in(:,:,n)';
        Ek_tmp=matrix(Ek_tmp,-1,1);
        Ek_out(:,n)=Ek_tmp;
    end
    Ek_out=cat(2,repmat([1:tot_k]',tot_ban,1),Ek_out);
    fid=mopen(file_name,'w');
    mfprintf(fid,'%4d'+strcat(repmat('   %10.6f',1,tot_proj))+'\n',Ek_out);
    mclose(fid);
endfunction
