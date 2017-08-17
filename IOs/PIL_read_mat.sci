// **** Purpose ****
// This is a read function of the standard output format of my library
// This format also consistent with my Fortran output format
// **** Variables ****
// [fid]: 1x1, integer
// <= your file ID
// [read_range]: 2x2, integer, default: read all data
// <= specify the range you want to read,ex:[3,10;7,12]
// [A]: NxM, integer, real or complex
// => you read data
// **** Version ****
// 02/07/2014 First Built
// 06/05/2014 auto read data type, no type input needed anymore
// 06/07/2014 add sparse and string matrix (separated by '#') format
// 04/26/2015 read row numbers 
// **** Comment ****
// 1. illustration
// Helpful to pass and receive data with Fortran. It makes Scilab and 
// Fortran well-integrated.

function [A]=PIL_read_mat(fid,read_range)
    // determin starting point and size
    while meof(fid)==0
        read_data=mgetl(fid,1);
        if length(grep(read_data,['PiLib Variable']))~=0 then
            mgetl(fid,1); 
            mgetstr(6,fid);
            order=mfscanf(1,fid,'%4d');
            mgetstr(8,fid);
            A_size(1)=mfscanf(1,fid,'%6d');
            mgetstr(1,fid);
            A_size(2)=mfscanf(1,fid,'%6d');
            mgetstr(8,fid);
            data_type=mfscanf(1,fid,'%s');
            if data_type=='STRING' then
                mgetl(fid,2);
            else
                mgetl(fid,4);
            end
            break;
        end
    end
    [lhs,rhs]=argn();
    if rhs==1 then
        read_range=[1,A_size(1);1,A_size(2)];
    end
    r_size=read_range(1,2)-read_range(1,1)+1;
    c_size=read_range(2,2)-read_range(2,1)+1;
    tot_column=A_size(2);
    A=zeros(r_size,c_size);

    mgetl(fid,read_range(1,1)-1);
    select data_type
    case 'INTEGER'
        for n=1:r_size
            mfscanf(1,fid,'%6d');
            mfscanf(read_range(2,1)-1,fid,'%d ')
            A(n,:)=(mfscanf(c_size,fid,'%d '))';
            if read_range(2,2)< tot_column then
                mgetl(fid,1);
            end
        end
    case 'REAL'
        for n=1:r_size
            mfscanf(1,fid,'%6d');
            mfscanf(read_range(2,1)-1,fid,'%f ')
            A(n,:)=(mfscanf(c_size,fid,'%f '))';
            if read_range(2,2)< tot_column then
                mgetl(fid,1);
            end
        end
        A=A*10^(order);
    case 'COMPLEX'
        A=zeros(r_size,c_size);
        B=zeros(r_size,c_size*2);
        for n=1:r_size
            mfscanf(1,fid,'%6d');
            mfscanf(2*(read_range(2,1)-1),fid,'%f ')
            B(n,:)=(mfscanf(2*c_size,fid,'%f '))';
            if read_range(2,2)< tot_column then
                mgetl(fid,1);
            end
        end

        for n=1:c_size
            A(:,n)=B(:,2*n-1)+%i*B(:,2*n)
        end
        A=A*10^(order);
    case 'SPARSE'
        if rhs~=1 then
            disp('Warning: PIL_read_mat, sparse format cannot assign read range!');
        end
        A=mfscanf(A_size(1),fid,'%f %f %f %f') 
        A(:,1:2)=round(A(:,1:2));
        A(:,3)=A(:,3)+%i*A(:,4);
        A=A(:,1:3);
        A(:,3)=A(:,3)*10^(order);
        A=sparse(PIL_sparse(A,'sparse','all'));
    case 'STRING'
        A=emptystr(A_size(1),A_size(2));
        for n=1:A_size(1)
            A(n,:)=PIL_str_split(mgetl(fid,1));
        end
        A=A(read_range(1,1):read_range(1,2),read_range(2,1):read_range(2,2));
    end
endfunction
